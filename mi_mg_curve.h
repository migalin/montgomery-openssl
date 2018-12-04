#include <stdio.h>
#include <openssl/bn.h>

struct mi_montgomery_curve {
    BIGNUM *A, *B, *C, *p, *q;
};

void mi_montgomery_curve_init(struct mi_montgomery_curve *curve);
void mi_montgomery_curve_mg_ladder(const struct mi_montgomery_curve *curve, struct mi_point *point, const BIGNUM *power);
void mi_montgomery_curve_free(struct mi_montgomery_curve *curve);
int mi_montgomery_curve_point_on_curve(const struct mi_montgomery_curve *curve, const struct mi_point *point);
