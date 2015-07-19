#include "xcint_blas.h"


extern "C" {
    extern void dgemm_(char         *ta,
                       char         *tb,
                       int          *m,
                       int          *n,
                       int          *k,
                       double       *alpha,
                       const double *a,
                       int          *lda,
                       const double *b,
                       int          *ldb,
                       double       *beta,
                       double       *c,
                       int          *ldc);
    extern void dsymm_(char         *si,
                       char         *up,
                       int          *m,
                       int          *n,
                       double       *alpha,
                       const double *a,
                       int          *lda,
                       const double *b,
                       int          *ldb,
                       double       *beta,
                       double       *c,
                       int          *ldc);
};


void xcint_dgemm(char         *ta,
                 char         *tb,
                 int          *m,
                 int          *n,
                 int          *k,
                 double       *alpha,
                 const double *a,
                 int          *lda,
                 const double *b,
                 int          *ldb,
                 double       *beta,
                 double       *c,
                 int          *ldc)
{
    dgemm_(ta,
           tb,
           m,
           n,
           k,
           alpha,
           a,
           lda,
           b,
           ldb,
           beta,
           c,
           ldc);
}


void xcint_dsymm(char         *si,
                 char         *up,
                 int          *m,
                 int          *n,
                 double       *alpha,
                 const double *a,
                 int          *lda,
                 const double *b,
                 int          *ldb,
                 double       *beta,
                 double       *c,
                 int          *ldc)
{
    dsymm_(si,
           up,
           m,
           n,
           alpha,
           a,
           lda,
           b,
           ldb,
           beta,
           c,
           ldc);
}
