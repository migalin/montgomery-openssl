/*
 * @author Sergey Migalin
 * @year   2018
 */

#include <stdio.h>        // IO
#include <openssl/bn.h>   // BIGNUM lib from openssl
#include "mi_params.h"    // Lab parameters
#include "mi_point.h"     // Point functions
#include "mi_mg_curve.h"  // Elliptic curve functions

int main ()
{
	printf("\n\e[36mOpenSSL (libcrypto) Montgomery elliptic curve lab\e[0m\n\n");

	BIGNUM *power = BN_new(),    						// power of point
		   *test_x = BN_new(),   						// var to test algorythm
		   *max_rand = NULL,     						// range of power randint
		   *one = NULL;          						// 1

	struct mi_point point = {NULL, NULL};
    struct mi_montgomery_curve curve = {NULL, NULL, NULL, NULL, NULL};

	BN_dec2bn(&one, "1");               				// make some "one" number
	mi_montgomery_curve_init(&curve);   				// initialize elliptic curve with parameters
	BN_dec2bn(&max_rand, MI_MAX_RAND);  				// initialize random range


	mi_point_std(&point);               				// initialize id-tc26-gost-3410-2012-256-paramSetA point
	printf("[INFO] X = %s,\n       Z = %s\n", BN_bn2dec(point.x), BN_bn2dec(point.z));


	BN_rand_range(power, max_rand);     				// make cryptorandom integer from 0 to max_rand -- power
	//BN_dec2bn(&power, "18791202894425416709884867221754140282839356181413286370677768873537813477604");
	printf("[INFO] Random power is %s\n", BN_bn2dec(power));

	
	mi_montgomery_curve_mg_ladder(&curve, &point, power);
	printf("\e[32m[ANSWER] X ** %s =\n         = %s\e[0m\n", BN_bn2dec(power), BN_bn2dec(point.x));
	printf("[INFO] X = %s,\n       Z = %s\n", BN_bn2dec(point.x), BN_bn2dec(point.z));


	printf("\e[32m[CHECK 1] Point is ");
	if (mi_montgomery_curve_point_on_curve(&curve, &point) != 1)
		printf("NOT ");
	printf("on curve\e[0m\n");


	BN_copy(test_x, point.x);
	BN_add(power, curve.q, one);
	mi_montgomery_curve_mg_ladder(&curve, &point, power);
	printf("\e[32m[CHECK 2] Point with Q+1 power is ");
	if (BN_cmp(point.x, test_x) != 0)
		printf("NOT ");
	printf("the same!\e[0m\n");
	printf("[INFO] X = %s,\n       Z = %s\n", BN_bn2dec(point.x), BN_bn2dec(point.z));


	printf("\e[32m[CHECK 3] Point with Q power is ");
	mi_montgomery_curve_mg_ladder(&curve, &point, curve.q);
	if (!BN_is_one(point.x))
		printf("NOT ");
	printf("neutral\e[0m\n");
	printf("[INFO] X = %s,\n       Z = %s\n", BN_bn2dec(point.x), BN_bn2dec(point.z));


	BN_free(power); BN_free(test_x);
	mi_point_free(&point);
	mi_montgomery_curve_free(&curve);


}
