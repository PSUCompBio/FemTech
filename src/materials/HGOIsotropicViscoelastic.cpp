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
 * \lambda = properties(2)
 * k_1     = properties(3)
 * k_2     = properties(4)
 * g_1     = properties(5)
 * t_1     = properties(6)
 * g_2     = properties(7)
 * t_2     = properties(8)
 * Implementation follows : Kaliske, M. and Rothert, H., 1997. Formulation and 
 * implementation of three-dimensional viscoelasticity at small and finite 
 * strains. Computational Mechanics, 19(3), pp.228-239. */

void HGOIsotropicViscoelastic(int e, int gp) {
  if (ndim == 2) {
    FILE_LOG_SINGLE(ERROR, "Plane Strain implementation yet to be done");
    exit(EXIT_FAILURE);
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
  const double lambda = localProperties[2];
  const double k1 = localProperties[3];
  const double k2 = localProperties[4];
  // Kappa = 1/3 for isotropic fiber distribution
  const double kappa = 1.0/3.0;
  const double g1 = localProperties[5];
  const double t1 = localProperties[6];
  const double g2 = localProperties[7];
  const double t2 = localProperties[8];

  // Compute the bulk modulus from \lambda and \mu
  const double K = lambda + 2.0*mu/3.0;
  double * const F_element_gp = &(F[index]);
  // Compute the hydrostatic diagonal contribution to cauchy stress tensor
  const double hydroDiag = 0.5*K*(J*J-1.0)/J; 
  // Compute the left elastic Cauchy-Green tensor B = F*F^T
  double matSize = ndim * ndim;
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

  // Update deviatoric part of S_0 for next time step
  for (int i = 0; i < matSize; ++i) {
    S0nLocal[i] = Sdev[i];
  }

#ifdef DEBUG
  FILE_LOG_SINGLE(DEBUGLOGIGNORE, "Element ID = %d, gp = %d", e, gp);
  for(int i=0;i<6;i++){
    int indexD = pk2ptr[e]+6*gp+i;
    FILE_LOG_SINGLE(DEBUGLOGIGNORE, "PK2[%d] = %10.6e", indexD, pk2[indexD]);
  }
#endif //DEBUG
}
