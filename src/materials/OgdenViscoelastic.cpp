#include "FemTech.h"
#include "blas.h"
#include "lapack.h"

#include <stdlib.h>

/* Implements the Ogden Model
 * Expected material properties :
 * \rho      = properties(0)
 * K         = properties(1)
 * N_{Ogden} = properties(2), maximum value of 3
 * \alpha_i  = properties(3+2*i) i = 0 to N_{Ogden}-1
 * \mu_i     = properties(4+2*i) i = 0 to N_{Ogden}-1
 * N_{Prony} = properties(3+N_{Ogden}*2), maximum value of 6
 * g_i       = properties(4+2*N_{Ogden}+2*j) j = 0 to N_{Prony}-1
 * \tau_i    = properties(5+2*N_{Ogden}+2*j) j = 0 to N_{Prony}-1
 * MAXMATPARAMS = 4+2*(N_{Ogden}+N_{Prony}) = 22
 * */

void Ogden(int e, int gp) {
  if (ndim == 2) {
    FILE_LOG_SINGLE(ERROR, "Plane Strain implementation yet to be done");
    TerminateFemTech(3);
  }
  // Assumes ndim == 3
  // Pointer to start of deformation gradient matrix for given element number
  // and Gauss point
  const int index = fptr[e] + ndim * ndim * gp;
  double J = detF[detFptr[e] + gp];
  // Get material property ID to read material properties
  const int pideIndex = pid[e]*MAXMATPARAMS;
  // Get location of array with material properties of the element
  const double* localProperties = &(properties[pideIndex]);
  // Read material properties
  const double K = localProperties[1];
  const int nTerm = static_cast<int>(localProperties[2]);

  double * const F_element_gp = &(F[index]);
  // Use temp storage mat1 for Finverse
  double *fInv = mat1;
  // Compute F inverse
  InverseF(e, gp, fInv);
  // Use temp storage 2 for storing Bmat
  double *Bmat = mat2;
  // Compute B = FF^T and compute its eigen values
  dgemm_(chn, chy, &ndim, &ndim, &ndim, &one, F_element_gp, &ndim,
        F_element_gp, &ndim, &zero, Bmat, &ndim);
  double matSize = ndim * ndim;
  // Compute eigen values and eigenvectors of B
  int nEigen;
  double cEigenValue[ndim];
  double cEigenVector[ndim*ndim];
  int dWorkN = ndim*30;
  int iWorkN = 15*ndim;
  double dWork[dWorkN];
  int iWork[iWorkN];
  int iSuppZ[2*ndim];
  dsyevr_(jobzV, rangeA, uploU, &ndim, Bmat, &ndim, &dStart, &dStart, \
      &iStart, &iStart, &eigenTol, &nEigen, cEigenValue, cEigenVector, \
      &ndim, iSuppZ, dWork, &dWorkN, iWork, &iWorkN, &info);
  // Compute sqrt of lambda and multiply by J^{-1/3}
  const double Jm13 = pow(J, -1.0/3.0);
  for (int i = 0; i < ndim; ++i) {
    cEigenValue[i] = sqrt(cEigenValue[i])*Jm13;
  }
  // Compute F^-1(cEigenVector)
  // Each column of basisVec points to F^{-1}*e where e is the eigen vector of B
  // matrix
  double *basisVec = mat3;
  dgemm_(chn, chn, &ndim, &ndim, &ndim, &one, fInv, &ndim,
           cEigenVector, &ndim, &zero, basisVec, &ndim);

  double *S = mat4;
  for (int i = 0; i < matSize; ++i) {
    S[i] = 0.0;
  }
  const double hydro = K*(J-1.0); 
  double eigenPower[ndim];
  for (int i = 0; i < ndim; ++i) {
    double preFactor = 0.0;
    for (int j = 0; j < nTerm; ++j) {
      double alpha = localProperties[3+j*2];
      double mu = localProperties[4+j*2];
      double eigenSum = 0.0;
      for (int k = 0; k < ndim; ++k) {
        const double t1 = pow(cEigenValue[k], alpha);
        eigenPower[k] = t1;
        eigenSum = eigenSum + t1;
      }
      eigenSum = eigenSum/(-3.0);
      preFactor = mu*(eigenPower[i]-eigenSum);
    }
    preFactor += hydro;
    dyadic(&basisVec[3*i], preFactor, S);
  }
  double *Cmat = mat2;
  dgemm_(chy, chn, &ndim, &ndim, &ndim, &one, F_element_gp, &ndim,
        F_element_gp, &ndim, &zero, Cmat, &ndim);
  // compute the double dot product C:S
  double SddC = 0.0;
  for (int i = 0; i < matSize; ++i) {
    SddC += Cmat[i]*S[i];
  }
  // Compute isochoric part : Sic = 1/3(C:S)inv(C)
  SddC = SddC/3.0;
  // Compute Inverse of C, inv(C) = inv(F)*inv(F)^T
  double* Sic = mat3; // Reuse mat3 for storing isochoric part of S
  dgemm_(chn, chy, &ndim, &ndim, &ndim, &SddC, fInv, &ndim,
        fInv, &ndim, &zero, Sic, &ndim);
  double *Sdev = mat1; // Reuse mat1 to store deviatoric part of S
  for (int i = 0; i < matSize; ++i) {
    Sdev[i] = S[i]-Sic[i];
  }
  // Access histroy dependence
  double *Hn_1Local = &(Hn_1[index]);
  double *Hn_2Local = &(Hn_2[index]);
  double *S0nLocal = &(S0n[index]);
  // Compute c1i = exp(-dt/t_i) and C2i = g_i (1-exp(-dt/t_i))/(dt/t_i)
  const double rt1 = dt/t1;
  const double rt2 = dt/t2;
  const double c11 = exp(-rt1);
  const double c12 = exp(-rt2);
  const double c21 = g1*(1-c11)/rt1;
  const double c22 = g2*(1-c12)/rt2;
  // Update H_j^{n+1} = c1j*H_j^n + c2j*[Dev S_0^{n+1} - S_0^n]
  // Update S^{n+1} = S_0^{n+1} + \sum_j H_j^{n+1}
  for (int i = 0; i < matSize; ++i) {
    Hn_1Local[i] = c11*Hn_1Local[i]+c21*(Sdev[i]-S0nLocal[i]);
    Hn_2Local[i] = c12*Hn_2Local[i]+c22*(Sdev[i]-S0nLocal[i]);
    S[i] = S[i] + Hn_1Local[i] + Hn_2Local[i];
  }

  // Get location of array to store PK2 values
  double * pk2Local = &(pk2[pk2ptr[e]+6*gp]);
	// 6 values saved per gauss point for 3d
	// in voigt notation, sigma11
  pk2Local[0] = S[0];
  // in voigt notation, sigma22
  pk2Local[1] = S[4];
  // in voigt notation, sigma33
  pk2Local[2] = S[8];
  // in voigt notation, sigma23
  pk2Local[3] = S[7];
  // in voigt notation, sigma13
  pk2Local[4] = S[6];
  // in voigt notation, sigma12
  pk2Local[5] = S[3];

#ifdef DEBUG
  FILE_LOG_SINGLE(DEBUGLOGIGNORE, "Element ID = %d, gp = %d", e, gp);
  for(int i=0;i<6;i++){
    int indexD = pk2ptr[e]+6*gp+i;
    FILE_LOG_SINGLE(DEBUGLOGIGNORE, "PK2[%d] = %10.6e", indexD, pk2[indexD]);
  }
#endif //DEBUG
}
