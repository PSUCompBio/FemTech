#include "FemTech.h"
#include "blas.h"

#include <stdlib.h>

/* Implements the three-dimensional Holzapfel-Gasser-Ogden Model with isotropic
 * fiber distribution (\kappa = \frac{1}{3}). Details of this implementation can
 * be found in HGOIsotropic function. Here on top of the HGO model, 
 * viscoelasticity is included for the shear compnents of the stress tensor.
 * Computes the second Piola-Kirchhoff stress tensor in Voigt notation. 
 * Input to the function is the eleemnt ID and the Gauss Quadrature point
 * number. Using the input the deformation gradient is loaded along with the
 * material properties of the element. Expected material properties :
 * \rho    = properties(0)
 * \mu     = properties(1)
 * K       = properties(2)
 * k_1     = properties(3)
 * k_2     = properties(4)
 * n       = properties(5)
 * g_i     = properties(6+2*i)
 * t_i     = properties(7+2*i)
 * Implementation follows : Kaliske, M. and Rothert, H., 1997. Formulation and 
 * implementation of three-dimensional viscoelasticity at small and finite 
 * strains. Computational Mechanics, 19(3), pp.228-239. */

void HGOIsotropicViscoelastic(int e, int gp) {
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
  const double mu = localProperties[1];
  const double K = localProperties[2];
  const double k1 = localProperties[3];
  const double k2 = localProperties[4];
  // Kappa = 1/3 for isotropic fiber distribution
  const double kappa = 1.0/3.0;
  const int nPronyLocal = int(localProperties[5]); // No of terms in prony series
  // TODO(Anil) convert local allocations to one time global
  double *gi, *ti;
  gi = (double*)malloc(nPronyLocal*sizeof(double));
  ti = (double*)malloc(nPronyLocal*sizeof(double));
  for (unsigned int i = 0; i < nPronyLocal; ++i) {
    gi[i] = localProperties[6+2*i];
    ti[i] = localProperties[7+2*i];
  }

  double * const F_element_gp = &(F[index]);
  // Compute the hydrostatic diagonal contribution to cauchy stress tensor
  const double hydroDiag = 0.5*K*(J*J-1.0)/J; 
  // Compute the left elastic Cauchy-Green tensor B = F*F^T
  const unsigned int matSize = ndim * ndim;
  double *Bmat = mat1;
  dgemm_(chn, chy, &ndim, &ndim, &ndim, &one, F_element_gp, &ndim,
        F_element_gp, &ndim, &zero, Bmat, &ndim);
  // Compute trace of B
  const double traceB = Bmat[0] + Bmat[4] + Bmat[8];
  // Calculate first invarient of the isochoric Cauchy-Green deformation tensor
  const double Jm23 = pow(J, -2.0/3.0);
  const double I1 = Jm23*traceB;
  // Compute E_\alpha
  const double E_alpha = kappa*(I1-3.0);
  double fiberPrefactor = 0.0;
  if (E_alpha > 0.0) {
    fiberPrefactor = 2.0*k1*exp(k2*E_alpha*E_alpha)*E_alpha*kappa;
  }
  const double totalPrefactor = (mu + fiberPrefactor)/J;
  // Compute deviatoric part of isochoric left elastic Cauchy-Green tensor
  // \overline{B}^d_{ij} = J^{-2/3}(B_{ij}-B_{kk} \delta_{ij}/3)
  const double traceBby3 = traceB/3.0;
  Bmat[0] = Bmat[0] - traceBby3;
  Bmat[4] = Bmat[4] - traceBby3;
  Bmat[8] = Bmat[8] - traceBby3;
  for (int i = 0; i < matSize; ++i) {
    Bmat[i] = Bmat[i]*Jm23*totalPrefactor;
  }
  // Bmat now stores \sigma_d
  // Storing sigma to Bmat
  Bmat[0] += hydroDiag;
  Bmat[4] += hydroDiag;
  Bmat[8] += hydroDiag;

  // Compute pk2 : S = J F^{-1} \sigma  F^{-T}
  double *fInv = mat2;
  double *STemp = mat3;
  double *S = mat4;
  InverseF(e, gp, fInv);
  // Compute F^{-1}*\sigma
  dgemm_(chn, chn, &ndim, &ndim, &ndim, &one, fInv, &ndim,
           Bmat, &ndim, &zero, STemp, &ndim);
  // Compute J F^{-1}*\sigma F^{-T}
  dgemm_(chn, chy, &ndim, &ndim, &ndim, &J, STemp, &ndim,
           fInv, &ndim, &zero, S, &ndim);
  // Compute the right elastic Cauchy-Green tensor C = F^T*F
  // Reuse temperory variable mat1 to compute C
  double *Cmat = mat1;
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
  double *S0nLocal = S0n[e];
  double *elemH = Hn[e];
  double *HnLocal = &(elemH[gp*matSize*nPronyLocal]);
  for (unsigned int i = 0; i < nPronyLocal; ++i) {
    double *HnI = &(HnLocal[i*matSize]);
    // Compute C1i = exp(-dt/t_i) and C2i = g_i (1-exp(-dt/t_i))/(dt/t_i)
    const double rt = dt/ti[i];
    const double c1i = exp(-rt);
    const double c2i = gi[i]*(1-c1i)/rt;
    // Update H_j^{n+1} = c1j*H_j^n + c2j*[Dev S_0^{n+1} - S_0^n]
    // Update S^{n+1} = S_0^{n+1} + \sum_j H_j^{n+1}
    for (int j = 0; j < matSize; ++j) {
      HnI[j] = c1i*HnI[j]+c2i*(Sdev[j]-S0nLocal[j]);
      S[j] = S[j] + HnI[j];
    }
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

  // Update deviatoric part of S_0 for next time step
  for (int i = 0; i < matSize; ++i) {
    S0nLocal[i] = Sdev[i];
  }
  free(gi);
  free(ti);

#ifdef DEBUG
  FILE_LOG_SINGLE(DEBUGLOGIGNORE, "Element ID = %d, gp = %d", e, gp);
  for(int i=0;i<6;i++){
    int indexD = pk2ptr[e]+6*gp+i;
    FILE_LOG_SINGLE(DEBUGLOGIGNORE, "PK2[%d] = %10.6e", indexD, pk2[indexD]);
  }
#endif //DEBUG
}
