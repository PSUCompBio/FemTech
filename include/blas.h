#ifndef INCLUDE_BLAS_H_
#define INCLUDE_BLAS_H_

// Variables for blas
static char* chy = (char*)"T";
static char* chn = (char*)"N";
static double one = 1.0;
static double zero = 0.0;
static int oneI = 1;

extern "C" {
//  extern double ddot_(int *n, double *dx, int *incx, double *dy, int *incy);
 extern void dgemm_(char *transa, char *transb, int *m, int *n, int *k, \
     double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, \
     double *c, int *ldc);
 extern void dgesv_(int *n, int *nrhs, double *a, int *lda, \
     int *ipiv, double *b, int *ldb, int *info);
 extern void dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
 extern void dgetrs_(char *trans, int *n, int *nrhs, double *a, int *lda, \
     int *ipiv, double *b, int *ldb, int *info);
 extern void dgemv_(char *trans, int *m, int *n, double *alpha, double *a, \
     int *lda, double *x, int *incx, double *beta, double *y, int *incy);
 extern void daxpy_(int *n, double *alpha, double *x, int *incx, double *y, int *incy);
}

#endif // INCLUDE_BLAS_H_
