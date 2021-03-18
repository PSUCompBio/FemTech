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
  const unsigned int matSize = ndim * ndim;
  // Compute eigen values and eigenvectors of B
  double cEigenValue[ndim];
  double dWork[workSize];
  dsyev_(jobzV, uploU, &ndim, Bmat, &ndim, cEigenValue, dWork, &workSize, &info);
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
           Bmat, &ndim, &zero, basisVec, &ndim);

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

  // TODO (NEED UPDATE from PK2 to Cauchy)
  // Get location of array to store Cauchy values
  double * sigma_nLocal = &(sigma_n[sigmaptr[e]+6*gp]);
  // 6 values saved per gauss point for 3d
  // in voigt notation, sigma11
  sigma_nLocal[0] = S[0];
  // in voigt notation, sigma22
  sigma_nLocal[1] = S[4];
  // in voigt notation, sigma33
  sigma_nLocal[2] = S[8];
  // in voigt notation, sigma23
  sigma_nLocal[3] = S[7];
  // in voigt notation, sigma13
  sigma_nLocal[4] = S[6];
  // in voigt notation, sigma12
  sigma_nLocal[5] = S[3];

#ifdef DEBUG
  FILE_LOG_SINGLE(DEBUGLOGIGNORE, "Element ID = %d, gp = %d", e, gp);
  for(int i=0;i<6;i++){
    int indexD = sigmaptr[e]+6*gp+i;
    FILE_LOG_SINGLE(DEBUGLOGIGNORE, "Sigma[%d] = %10.6e", indexD, sigma_n[indexD]);
  }
#endif //DEBUG
}
