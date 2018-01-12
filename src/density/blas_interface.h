#pragma once

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
                int *ldc);

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
                int *ldc);
