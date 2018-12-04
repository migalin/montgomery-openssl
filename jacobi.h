#include <openssl/bn.h>

#define BN_lsw(n) (((n)->top == 0) ? (BN_ULONG) 0 : (n)->d[0])

int bn_jacobi(const BIGNUM *a, const BIGNUM *b, BN_CTX *ctx) {
  // In 'tab', only odd-indexed entries are relevant:
  // For any odd BIGNUM n,
  //     tab[BN_lsw(n) & 7]
  // is $(-1)^{(n^2-1)/8}$ (using TeX notation).
  // Note that the sign of n does not matter.
  static const int tab[8] = {0, 1, 0, -1, 0, -1, 0, 1};

  int ret = -2;
  BN_CTX_start(ctx);
  BIGNUM *A = BN_CTX_get(ctx);
  BIGNUM *B = BN_CTX_get(ctx);
  if (B == NULL) {
    goto end;
  }

  if (!BN_copy(A, a) ||
      !BN_copy(B, b)) {
    goto end;
  }
  ret = 1;

  while (1) {
    // B is positive and odd
    if (BN_is_zero(A)) {
      ret = BN_is_one(B) ? ret : 0;
      goto end;
    }
    // now A is non-zero
    int i = 0;
    while (!BN_is_bit_set(A, i)) {
      i++;
    }
    if (!BN_rshift(A, A, i)) {
      ret = -2;
      goto end;
    }
    if (i & 1) {
      // i is odd
      // multiply 'ret' by  $(-1)^{(B^2-1)/8}$
      ret = ret * tab[BN_lsw(B) & 7];
    }
    // multiply 'ret' by  $(-1)^{(A-1)(B-1)/4}$
    if ((A->neg ? ~BN_lsw(A) : BN_lsw(A)) & BN_lsw(B) & 2) {
      ret = -ret;
    }
    // (A, B) := (B mod |A|, |A|)
    if (!BN_nnmod(B, B, A, ctx)) {
      ret = -2;
      goto end;
    }
    BIGNUM *tmp = A;
    A = B;
    B = tmp;
    tmp->neg = 0;
  }
end:
  BN_CTX_end(ctx);
  return ret;
}
