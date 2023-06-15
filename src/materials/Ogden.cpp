#include "FemTech.h"
#include "blas.h"
#include "lapack.h"

#include <stdlib.h>
#include <math.h>

/* Implements the Ogden Model
 * Expected material properties :
 * \rho      = properties(0)
 * K         = properties(1)
 * N_{Ogden} = properties(2), maximum value of 3
 * \alpha_i  = properties(3+2*i) i = 0 to N_{Ogden}-1
 * \mu_i     = properties(4+2*i) i = 0 to N_{Ogden}-1
 * MAXMATPARAMS = 4+2*N_{Ogden} = 10
 * */

/*
 * compute Cauchy Stress \sigma
 * \sigma = \Sum_{i=0}^{3} ((\Sum_{j=0}^{N} mu_j [\overline{\lambda}_i^{\alpha_j}
 * - (\overline{\lambda}_1^{\alpha_j}+\overline{\lambda}_2^{\alpha_j}+
 *   \overline{\lambda}_3^{\alpha_j})/3.0])) (n_i \dyadic n_i) + K*(J-1)/J
 */

void Ogden(int e, int gp) {
  if (ndim == 2) {
    FILE_LOG_SINGLE(ERROR, "Plane Strain implementation yet to be done");
    TerminateFemTech(3);
  }
  // Assumes ndim == 3
  // Get material property ID to read material properties
  const int pideIndex = pid[e]*MAXMATPARAMS;
  // Get location of array with material properties of the element
  const double* localProperties = &(properties[pideIndex]);
  // Read material properties
  const double K = localProperties[1];
  const int nTerm = static_cast<int>(localProperties[2]);

  // b = F F^T
  double *H = mat1;
  ComputeH(e, gp, H);
  // Compute and store F = H + I
  const unsigned int index = fptr[e] + ndim2 * gp;
  double * const F_element_gp = &(F[index]);
  for (unsigned int i = 0; i < ndim2; ++i) {
    F_element_gp[i] = H[i];
  }
  F_element_gp[0] = F_element_gp[0] + 1.0;
  F_element_gp[4] = F_element_gp[4] + 1.0;
  F_element_gp[8] = F_element_gp[8] + 1.0;
  double *b = mat2;
  // Compute F F^T and store to b 
  dgemm_(chn, chy, &ndim, &ndim, &ndim, &one, F_element_gp, &ndim,
          F_element_gp, &ndim, &zero, b, &ndim);
  // Compute eigen values and eigenvectors of B
  //double cEigenValue[ndim];
  double* cEigenValue = new double[ndim];
  //double dWork[workSize];
  double* dWork = new double[workSize];
  dsyev_(jobzV, uploU, &ndim, b, &ndim, cEigenValue, dWork, &workSize, &info);
  // b now contains the eigen vectors of b

  // Correct Eigen Values for numerical inaccuracy
  long double delat01 = fabsl(cEigenValue[0] - cEigenValue[1]);
  long double delat12 = fabsl(cEigenValue[1] - cEigenValue[2]);
  long double delat02 = fabsl(cEigenValue[0] - cEigenValue[2]);

  const double eigenTol = 1e-12;
  if (delat01 < eigenTol) cEigenValue[1] = cEigenValue[0];
  if (delat02 < eigenTol) cEigenValue[2] = cEigenValue[0];
  if (delat12 < eigenTol) cEigenValue[2] = cEigenValue[1];

  const double J = det3x3Matrix(F_element_gp);
  const double Jm13 = pow(J, -1.0/3.0);
  const double invJ = 1.0/J;

  //double eigenPower[ndim*nTerm];
  double* eigenPower = new double[ndim * nTerm];
  //double eigenPowerSum[nTerm];
  double* eigenPowerSum = new double[nTerm];
  for (int j = 0; j < nTerm; ++j) {
    eigenPowerSum[j] = 0.0;
  }
  for (int i = 0; i < ndim; ++i) {
    const double lambdaBar = Jm13*sqrt(cEigenValue[i]);
    for (int j = 0; j < nTerm; ++j) {
      const double alpha = localProperties[3+j*2];
      const double lambdaBarAlpha = pow(lambdaBar, alpha);
      eigenPowerSum[j] += lambdaBarAlpha;
      eigenPower[i*nTerm+j] = lambdaBarAlpha;
    }
  }
  for (int j = 0; j < nTerm; ++j) {
    eigenPowerSum[j] /= 3.0;
  }

  double *sigma_e = mat3;
  for (int i = 0; i < ndim2; ++i) {
    sigma_e[i] = 0.0;
  }
  for (int i = 0; i < ndim; ++i) {
    double preFactor = 0.0;
    for (int j = 0; j < nTerm; ++j) {
      double mu = localProperties[4+j*2];
      preFactor += mu*(eigenPower[i*nTerm+j]-eigenPowerSum[j]);
    }
    preFactor *= invJ;
    dyadic(&b[3*i], preFactor, sigma_e);
  }
  const double volum = K*invJ*(J-1.0); 
  // const double volum = K*(J-1.0); 
  // Add the volumentric component
  sigma_e[0] += volum;
  sigma_e[4] += volum;
  sigma_e[8] += volum;

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

  delete[]cEigenValue;/*Drupal*/
  delete[]dWork;/*Drupal*/
  delete[]eigenPower;/*Drupal*/
  delete[]eigenPowerSum;/*Drupal*/

#ifdef DEBUG
  FILE_LOG_SINGLE(DEBUGLOGIGNORE, "Element ID = %d, gp = %d", e, gp);
  for(int i=0;i<6;i++){
    int indexD = sigmaptr[e]+6*gp+i;
    FILE_LOG_SINGLE(DEBUGLOGIGNORE, "Sigma[%d] = %10.6e", indexD, sigma_n[indexD]);
  }
#endif //DEBUG
}
