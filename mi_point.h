#include <openssl/bn.h>

struct mi_point {
	BIGNUM *x, *z;
};

void mi_point_neutral(struct mi_point *point);
void mi_point_std(struct mi_point *point);
void mi_point_add(struct mi_point *q, const struct mi_point *r, const struct mi_point *p1, const BIGNUM *p);
void mi_point_double(struct mi_point *point, const BIGNUM *c, const BIGNUM *p);
void mi_point_free(struct mi_point *point);

