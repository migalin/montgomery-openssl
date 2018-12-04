#include <openssl/bn.h>
#include "mi_params.h"

struct mi_point {
	BIGNUM *x, 
           *z;
};

/*
 * Initialize id-tc26-gost-3410-2012-256-paramSetA point
 * Make tranformation from Twisted Edwards to Montgomery parameters
 * 
 * Montgomery: b*y^2=x^3+a*x^2+x mod p
 * Twisted edrawds: eu^2 +v^2 =1+du^2v^2 mod p
 * 
 * Transform: x = (1 + v) / (1 - v)
 */
void mi_point_std(struct mi_point *point){

	BIGNUM *p = NULL,
		   *v = NULL,
	       *g = NULL,
	       *a = BN_new(),
	       *f = BN_new();

	BN_CTX *tmp = BN_CTX_new();           // tmp var for openssl

	BN_hex2bn(&p, MI_MG_MODULE);          // module
	BN_hex2bn(&v, MI_V_PARAMETER);        // v parameter of Twisted Edwards
	BN_dec2bn(&g, "1");                   // some "one"
	BN_dec2bn(&point->z, "1");            // another "one"

	BN_sub(f, g, v);                      // f = 1 - v
	BN_add(g, g, v);                      // g = 1 + v
	BN_mod_inverse(a, f, p, tmp);         // find invert element to (1 - v) in GF(p)
	if (!point->x) point->x = BN_new();   // if point was not initialized earlier
	BN_mod_mul(point->x, a, g, p, tmp);   // x = a / g mod p
	
	BN_free(a); BN_free(p); BN_free(v); BN_free(g); BN_free(f); BN_CTX_free(tmp);
}

/*
 * Initialize neutral point
 * On pojective coordinates 0 = (0,1,0)
 * Let x=1, z=0
 */
void mi_point_neutral(struct mi_point *point){
	BN_dec2bn(&point->x, "1");	
	BN_dec2bn(&point->z, "0");
}


/*
 * Add two points q and r
 * Transform: x_res = (((q.x - q.z) * (r.x + r.z) + (q.x + q.z) * (r.x - r.z))^2) * z1
 *            z_res = (((q.x - q.z) * (r.x + r.z) - (q.x + q.z) * (r.x - r.z))^2) * x1
 * Source: Daniel J. Bernstein and Tanja Lange "Montgomery curves and the Montgomery ladder" p.4
 */
void mi_point_add(struct mi_point *q, const struct mi_point *r, const struct mi_point *p1, const BIGNUM *p){
	BIGNUM *s_q  = BN_new(), 
		   *d_q = BN_new(), 
		   *s_r  = BN_new(),
		   *d_r = BN_new(),
		   *pow_2  = NULL;

    BN_CTX *tmp = BN_CTX_new();            // tmp var
	BN_dec2bn(&pow_2, "2");                // 2 number -- power

	BN_sub(d_q, q->x, q->z);               // d_q = q.x - q.z
    BN_add(s_q, q->x, q->z);               // s_q = q.x + q.z
	BN_add(s_r, r->x, r->z);               // s_r = r.x + r.z
	BN_sub(d_r, r->x, r->z);               // d_r = r.x - r.z

	BN_mul(d_q, d_q, s_r, tmp);            // d_q = (q.x - q.z) * (r.x + r.z)
	BN_mul(s_q, s_q, d_r, tmp);            // s_q = (q.x + q.z) * (r.x - r.z)

	BN_add(q->x, d_q, s_q);                // q.x = (q.x - q.z) * (r.x + r.z) + (q.x + q.z) * (r.x - r.z)
	BN_exp(q->x, q->x, pow_2, tmp);        // q.x = ((q.x - q.z) * (r.x + r.z) + (q.x + q.z) * (r.x - r.z))^2
	BN_mod_mul(q->x, q->x, p1->z, p, tmp); // q.x = ((q.x - q.z) * (r.x + r.z) + (q.x + q.z) * (r.x - r.z))^2 * p1.z mod p

	BN_sub(q->z, d_q, s_q);                // q.z = ((q.x - q.z) * (r.x + r.z)) - (q.x + q.z) * (r.x - r.z)
	BN_exp(q->z, q->z, pow_2, tmp);        // q.z = ((q.x - q.z) * (r.x + r.z)) - (q.x + q.z) * (r.x - r.z))^2
	BN_mod_mul(q->z, q->z, p1->x, p, tmp); // q.z = ((q.x - q.z) * (r.x + r.z)) - (q.x + q.z) * (r.x - r.z))^2 * p1.x mod p

    BN_free(s_q); BN_free(d_q); BN_free(s_r); BN_free(d_r); BN_free(pow_2);
	BN_CTX_free(tmp);
}


/*
 * Multiple point by 2
 * Transform: x_res = (point.x + point.z)^2 * (point.x - point.z)^2
 *            z_res = ((point.x + point.z)^2 - (point.x - point.z)^2) ((a - 2) / 4) * ((point.x + point.z)^2 - 
 *                  (point.x + point.z)^2) + (point.x + point.z)^2)
 * Source: Daniel J. Bernstein and Tanja Lange "Montgomery curves and the Montgomery ladder" p.4
 */
void mi_point_double(struct mi_point *point, const BIGNUM *c, const BIGNUM *p){

	BIGNUM *s_p   = BN_new(), 
		   *d_p  = BN_new(),
		   *mul_res = BN_new(),
		   *pow_2   = NULL;

	BN_CTX *tmp = BN_CTX_new();                    // tmp var
	BN_dec2bn(&pow_2, "2");                        // 2 number -- power

	BN_add(s_p, point->x, point->z);               // s_p = x + z
	BN_exp(s_p, s_p, pow_2, tmp);                  // s_p = (x + z)^2

	BN_sub(d_p, point->x, point->z);               // d_p = x - z
	BN_exp(d_p, d_p, pow_2, tmp);                  // d_p = (x - z)^2
	BN_mod_mul(point->x, s_p, d_p, p, tmp);        // point.x = (x + z)^2 * (x - z)^2 mod p
	BN_sub(d_p, s_p, d_p);                         // d_p = (x + z)^2 - (x - z)^2
	BN_mul(mul_res, c, d_p, tmp);                  // mul_res = c * d_p
	BN_add(s_p, s_p, mul_res);                     // s_p = (x + z)^2 + c * (x - z)^2
	BN_mod_mul(point->z, d_p, s_p, p, tmp);        // point.z = (x - z)^2 * (x + z)^2 + c * (x - z)^2 mod p

	BN_free(s_p); BN_free(d_p); BN_free(mul_res); BN_free(pow_2);
	BN_CTX_free(tmp);
}


/*
 * Clears memory
 */
void mi_point_free(struct mi_point *point){
	BN_free(point->x);
	BN_free(point->z);
}