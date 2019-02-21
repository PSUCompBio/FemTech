#ifndef INCLUDE_BLAS_H_
#define INCLUDE_BLAS_H_

// Variables for blas
static char* chy = (char*)"T";
static char* chn = (char*)"N";
static double one = 1.0;
static double zero = 0.0;

extern "C" {
//  extern double ddot_(int *n, double *dx, int *incx, double *dy, int *incy);
 extern void dgemm_(char *transa, char *transb, int *m, int *n, int *k, \
     double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, \
     double *c, int *ldc);
 extern void dgesv_(int *n, int *nrhs, double *a, int *lda, \
     int *ipiv, double *b, int *ldb, int *info);
}

#endif // INCLUDE_BLAS_H_
