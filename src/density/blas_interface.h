#pragma once

void wrap_dgemm(char ta,
                char tb,
                int m,
                int n,
                int k,
                double alpha,
                const double *a,
                int da,
                const double *b,
                int ldb,
                double beta,
                double *c,
                int ldc);

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
                int ldc);

extern "C"
{
    void dgemm_(const char *ta,
                const char *tb,
                const int *m,
                const int *n,
                const int *k,
                const double *alpha,
                const double *a,
                const int *lda,
                const double *b,
                const int *ldb,
                const double *beta,
                double *c,
                const int *ldc);
    void dsymm_(const char *si,
                const char *up,
                const int *m,
                const int *n,
                const double *alpha,
                const double *a,
                const int *lda,
                const double *b,
                const int *ldb,
                const double *beta,
                const double *c,
                int *ldc);
};
