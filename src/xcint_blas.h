#ifndef XCINT_BLAS_H_INCLUDED
#define XCINT_BLAS_H_INCLUDED

void xcint_dgemm(char *ta, char *tb, int *m, int *n, int *k, double *alpha,
                 const double *a, int *lda, const double *b, int *ldb,
                 double *beta, double *c, int *ldc);

void xcint_dsymm(char *si, char *up, int *m, int *n, double *alpha,
                 const double *a, int *lda, const double *b, int *ldb,
                 double *beta, double *c, int *ldc);

#endif /* XCINT_BLAS_H_INCLUDED */
