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

  // Get material property ID to read material properties
  const int pideIndex = pid[e]*MAXMATPARAMS;
  // Get location of array with material properties of the element
  const double* localProperties = &(properties[pideIndex]);
  // Read material properties
  const double K = localProperties[1];
  const int nTerm = static_cast<int>(localProperties[2]);
  const unsigned int matSize = ndim * ndim;
  double *H = mat1;
  ComputeH(e, gp, H);
  // Compute and store F = H + I
  double * const F_element_gp = &(F[index]);
  for (unsigned int i = 0; i < ndim2; ++i) {
      F_element_gp[i] = H[i];
      }
  F_element_gp[0] = F_element_gp[0] + 1.0;
  F_element_gp[4] = F_element_gp[4] + 1.0;
  F_element_gp[8] = F_element_gp[8] + 1.0;
  // Use temp storage mat2 for Finverse
  double *fInv = mat2;
  // Compute F inverse
  double J = inverse3x3Matrix(F_element_gp, fInv);
  //InverseF(e, gp, fInv);
  // Use temp storage 2 for storing Bmat
  double *Bmat = mat3;
  // Compute B = FF^T and compute its eigen values
  dgemm_(chn, chy, &ndim, &ndim, &ndim, &one, F_element_gp, &ndim,
        F_element_gp, &ndim, &zero, Bmat, &ndim);

   double invJ = 1.0/J;
  // Compute eigen values and eigenvectors of B
  int nEigen;
  double cEigenValue[ndim];
  double dWork[workSize];
  
    // Correct Eigen Values for numerical inaccuracy
  double delat01 = fabs(cEigenValue[0] - cEigenValue[1]);
  double delat12 = fabs(cEigenValue[1] - cEigenValue[2]);
  double delat02 = fabs(cEigenValue[0] - cEigenValue[2]);

  const double eigenTol = 1e-12;
  if (delat01 < eigenTol) cEigenValue[1] = cEigenValue[0];
  if (delat02 < eigenTol) cEigenValue[2] = cEigenValue[0];
  if (delat12 < eigenTol) cEigenValue[2] = cEigenValue[1];

  dsyev_(jobzV, uploU, &ndim, Bmat, &ndim, cEigenValue, dWork, &workSize, &info);
  // Compute sqrt of lambda and multiply by J^{-1/3}
  const double Jm13 = pow(J, -1.0/3.0);
  for (int i = 0; i < ndim; ++i) {
    cEigenValue[i] = sqrt(cEigenValue[i])*Jm13;
  }
  // Compute F^-1(cEigenVector)
  // Each column of basisVec points to F^{-1}*e where e is the eigen vector of B
  // matrix
  //double *basisVec = mat3;
  //dgemm_(chn, chn, &ndim, &ndim, &ndim, &one, fInv, &ndim,
    //       Bmat, &ndim, &zero, basisVec, &ndim);

  double *sigma_e = mat1; //check
  for(int i = 0; i<matSize; i++){
    sigma_e[i] = 0.0;
  }

  const double hydro = K*(J-1.0);
  double eigenPower[ndim];
  for (int i = 0; i < ndim; ++i) {
    double preFactor = hydro;
    for (int j = 0; j < nTerm; ++j) {
      double alpha = localProperties[3+j*2];
      double mu = localProperties[4+j*2];
      preFactor += (mu*invJ*(pow(cEigenValue[i], alpha)-((pow(cEigenValue[0], alpha)+pow(cEigenValue[1], alpha)+pow(cEigenValue[2], alpha))/3.0)));
    }
    dyadic(&Bmat[3*i], preFactor, sigma_e);
  }

  // Get location of array to store Cauchy values
  double * sigma_nLocal = &(sigma_n[sigmaptr[e]+6*gp]);
  // 6 values saved per gauss point for 3d
  // in voigt notation, sigma11
  sigma_nLocal[0] = sigma_e[0];
  // in voigt notation, sigma22
  sigma_nLocal[1] = sigma_e[4];
  // in voigt notation, sigma33
  sigma_nLocal[2] = sigma_e[8];
  // in voigt notation, sigma23
  sigma_nLocal[3] = sigma_e[7];
  // in voigt notation, sigma13
  sigma_nLocal[4] = sigma_e[6];
  // in voigt notation, sigma12
  sigma_nLocal[5] = sigma_e[3];

#ifdef DEBUG
  FILE_LOG_SINGLE(DEBUGLOGIGNORE, "Element ID = %d, gp = %d", e, gp);
  for(int i=0;i<6;i++){
    int indexD = sigmaptr[e]+6*gp+i;
    FILE_LOG_SINGLE(DEBUGLOGIGNORE, "Sigma[%d] = %10.6e", indexD, sigma_n[indexD]);
  }
#endif //DEBUG
}
