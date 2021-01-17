#include "FemTech.h"
#include "blas.h"
#include "lapack.h"

#include <stdlib.h>

/* Implements the Ogden Model
 * Expected material properties :
 * \rho      = properties(0)
 * K         = properties(1)
 * G_0       = properties(2)
 * G_\infity = properties(3)
 * \beta     = properties(4)
 * */

void Viscoelastic(int e, int gp) {
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
  const double mu0 = localProperties[3];
  const double twomu0 = 2.0*mu0;
  const double gamma = (localProperties[2]/mu0-1); // TODO : Precompute
  const double beta = localProperties[4];
  // const double lambda = K - 2.0*mu0/3.0;

  // Compute Green-Lagrange Tensor: eps = (1/2)*(F^T + F) - I
  const unsigned int matSize = ndim * ndim;
  // Use temp storage 1 for storing eps matrix
  double *eps = mat1;
  // Get element F matrix
  double *F_element_gp = &(F[index]);
  for (int i = 0; i < ndim; ++i) {
    for (int j = 0; j < ndim; ++j) {
      const int indexL = j + i*ndim;
      eps[indexL] = 0.5*(F_element_gp[indexL]+F_element_gp[i+j*ndim]);
    }
  }
  eps[0] -= 1.0;
  eps[4] -= 1.0;
  eps[8] -= 1.0;

  // Compute trace of eps matrix
  const double traceEps = eps[0]+eps[4]+eps[8];
  // Precompute 2*mu0*tr(eps)/3.0
  const double const_1 = 2.0*mu0*traceEps/3.0;
  // Precompute K*tr(eps)
  const double const_2 = K*traceEps;

  // Compute Cauchy Stress
  // sigma = lambda*tr(eps)*I+2*mu0*eps
  // Deviatoric part sigmaDev = 2*mu0*(eps-tr(eps)*I/3)
  // Isochoric part sigmaIso = K*tr(eps)*I
  // sigma = sigmaDev + sigmaIso
  // Use temp storage 2 for storing sigmadev matrix
  double *sigmaDev = mat2;
  for (int i = 0; i < matSize; ++i) {
    sigmaDev[i] = twomu0 * eps[i];
  }
  sigmaDev[0] -= const_1;
  sigmaDev[4] -= const_1;
  sigmaDev[8] -= const_1;

  // Reuse temp storage 1 for storing S matrix
  double *sigma = mat1;
  for (int i = 0; i < matSize; ++i) {
    sigma[i] = sigmaDev[i];
  }
  sigma[0] += const_2;
  sigma[4] += const_2;
  sigma[8] += const_2;

  // Access histroy dependence
  // Restore previous time step Deviatoric part of sigma
  double *elemS = S0n[e];
  double *S0nLocal = &(elemS[gp*matSize]);
  double *elemH = Hn[e];
  double *HnLocal = &(elemH[gp*matSize]);
  // Compute C1i = exp(-dt/t_i) and C2i = g_i (1-exp(-dt/t_i))/(dt/t_i)
  const double rt = dt*beta;
  const double c1 = exp(-rt);
  double c2;
  if (fabs(c1-1) < 1e-12) {
    c2 = 0.0;
  } else {
    c2 = gamma*(1-c1)/rt;
  }
  // Update H^{n+1} = c1*H_j^n + c2*[Dev S_0^{n+1} - Dev S_0^n]
  // Update sigma^{n+1} = sigma_0^{n+1} + H^{n+1}
  for (int j = 0; j < matSize; ++j) {
    // Along with computation, previous value of H is also updated in this
    // step
    HnLocal[j] = c1*HnLocal[j]+c2*(sigmaDev[j]-S0nLocal[j]);
    sigma[j] = sigma[j] + HnLocal[j];
  }
  // FILE_LOGMatrix_SINGLE(WARNING, HnLocal, 3, 3, "Gauss Point : %d, Hn", gp);
  FILE_LOGMatrix_SINGLE(DEBUGLOGIGNORE, sigma, 3, 3, "Gauss Point : %d, sigma", gp);
  // FILE_LOGMatrix(WARNING, sigma, ndim, ndim, "sigma Visco mat\n");

  // Compute pk2 : S = F^{-1} sigma 
  double *fInv = mat3;
  double *S = mat4;
  InverseF(e, gp, fInv);
  // Compute F^{-1}*sigma
  dgemm_(chn, chn, &ndim, &ndim, &ndim, &one, fInv, &ndim,
          sigma, &ndim, &zero, S, &ndim);

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

  // Update deviatoric part of S_0 for next time step
  for (int i = 0; i < matSize; ++i) {
    S0nLocal[i] = sigmaDev[i];
  }

#ifdef DEBUG
  FILE_LOG_SINGLE(DEBUGLOGIGNORE, "Element ID = %d, gp = %d", e, gp);
  for(int i=0;i<6;i++){
    int indexD = pk2ptr[e]+6*gp+i;
    FILE_LOG_SINGLE(DEBUGLOGIGNORE, "PK2[%d] = %10.6e", indexD, pk2[indexD]);
  }
#endif //DEBUG
}
