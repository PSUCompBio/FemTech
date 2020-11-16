#ifndef INCLUDE_LAPACK_H_
#define INCLUDE_LAPACK_H_

// Variables for blas
static const char* jobzV = (char*)"V";
static const char* rangeA = (char*)"A";
static const char* uploU = (char*)"U";
static const double dStart = 0;
static const int iStart = 0;
static const double eigenTol = 1e-12;
static int info;

extern "C" {
 /* DSYEVR prototype */
 extern void dsyevr_(const char* const jobz, const char* const range, 
     const char* const uplo, const int* const n, double* const a, 
     const int* const lda, const double* const vl, const double* const vu, 
     const int* const il, const int* const iu, const double* const abstol,
     int* const m, double* const w, double* const z, const int* const ldz, 
     int* const isuppz, double* const work, const int* const lwork, 
     int* const iwork, const int* const liwork, int* const info);
 extern void dger_(const int* const M, const int* const N,
     const double* const alpha, const double* const X, const int* const incX,
     const double* const Y, const int* const incY, double* const A, 
     const int* const lda);
}

#endif // INCLUDE_LAPACK_H_
