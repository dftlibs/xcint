#ifdef HAVE_MKL_BLAS
#include "mkl.h"
#endif

#include "blas_interface.h"

#ifndef HAVE_MKL_BLAS
extern "C" {
extern void dgemm_(char *ta,
                   char *tb,
                   int *m,
                   int *n,
                   int *k,
                   double *alpha,
                   const double *a,
                   int *lda,
                   const double *b,
                   int *ldb,
                   double *beta,
                   double *c,
                   int *ldc);
extern void dsymm_(char *si,
                   char *up,
                   int *m,
                   int *n,
                   double *alpha,
                   const double *a,
                   int *lda,
                   const double *b,
                   int *ldb,
                   double *beta,
                   double *c,
                   int *ldc);
};
#endif

void wrap_dgemm(char *ta,
                char *tb,
                int *m,
                int *n,
                int *k,
                double *alpha,
                const double *a,
                int *lda,
                const double *b,
                int *ldb,
                double *beta,
                double *c,
                int *ldc)
{
#ifdef HAVE_MKL_BLAS
    CBLAS_TRANSPOSE cta = (*ta == 't') ? CblasNoTrans : CblasTrans;
    CBLAS_TRANSPOSE ctb = (*tb == 't') ? CblasNoTrans : CblasTrans;
    cblas_dgemm(CblasRowMajor,
                cta,
                ctb,
                *m,
                *n,
                *k,
                *alpha,
                a,
                *lda,
                b,
                *ldb,
                *beta,
                c,
                *ldc);
#else
    dgemm_(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#endif
}

void wrap_dsymm(char *si,
                char *up,
                int *m,
                int *n,
                double *alpha,
                const double *a,
                int *lda,
                const double *b,
                int *ldb,
                double *beta,
                double *c,
                int *ldc)
{
#ifdef HAVE_MKL_BLAS
    CBLAS_SIDE cs = (*si == 'r') ? CblasLeft : CblasRight;
    CBLAS_UPLO cu = (*up == 'u') ? CblasLower : CblasUpper;
    cblas_dsymm(CblasRowMajor,
                cs,
                cu,
                *n,
                *m,
                *alpha,
                a,
                *lda,
                b,
                *ldb,
                *beta,
                c,
                *ldc);
#else
    dsymm_(si, up, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
#endif
}
