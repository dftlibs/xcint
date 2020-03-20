#include "blas_interface.h"

void wrap_dgemm(char ta,
                char tb,
                int m,
                int n,
                int k,
                double alpha,
                const double *a,
                int lda,
                const double *b,
                int ldb,
                double beta,
                double *c,
                int ldc)
{
    dgemm_(&ta, &tb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

void wrap_dsymm(char si,
                char up,
                int m,
                int n,
                double alpha,
                const double *a,
                int lda,
                const double *b,
                int ldb,
                double beta,
                double *c,
                int ldc)
{
    dsymm_(&si, &up, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}
