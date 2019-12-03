#include "FemTech.h"
#include "blas.h"

#include <stdlib.h>

/* Implements the three-dimensional Holzapfel-Gasser-Ogden Model with isotropic
 * fiber distribution (\kappa = \frac{1}{3}). Only one set of fiber is assumed.
 * Computes the second Piola-Kirchhoff stress tensor in Voigt notation. 
 * Input to the function is the eleemnt ID and the Gauss Quadrature point
 * number. Using the input the deformation gradient is loaded along with the
 * material properties of the element. Expected material properties :
 * \rho    = properties(0)
 * \mu     = properties(1)
 * \lambda = properties(2)
 * k_1     = properties(3)
 * k_2     = properties(4)
 * Follows the formulation in Abaqus 6.13 and uses the formulae from 
 * Micromechanics of diffuse axonal injury: influence of axonal orientation and anisotropy
 * Cloots et al. Biomech Model Mechanobiol (2011)[Ref 1] */

void HGOIsotropic(int e, int gp) {
  if (ndim == 2) {
    printf("Plane Strain implementation yet to be done \n");
    exit(EXIT_FAILURE);
  }
  // Assumes ndim == 3
  // Pointer to start of deformation gradient matrix for given element number
  // and Gauss point
  const int index = fptr[e] + ndim * ndim * gp;
  const double J = detF[detFptr[e] + gp];
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

  // Compute the bulk modulus from \lambda and \mu
  const double K = lambda + 2.0*mu/3.0;
  double * const F_element_gp = &(F[index]);
  // Compute the hydrostatic diagonal contribution to cauchy stress tensor
  const double hydroDiag = 0.5*K*(J*J-1.0)/J; 
  // Compute the left elastic Cauchy-Green tensor B = F*F^T
  double matSize = ndim * ndim;
  double *Bmat = (double *)malloc(matSize * sizeof(double));
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
    Bmat[i] = Bmat[i]*Jm23;
  }

  // Get location of array to store PK2 values
  double * pk2Local = &(pk2[pk2ptr[e]+6*gp]);

	// in voigt notation, sigma11
	pk2Local[0] = hydroDiag + fiberPrefactor*Bmat[0];
	// in voigt notation, sigma22
	pk2Local[1] = hydroDiag + fiberPrefactor*Bmat[4];
	// in voigt notation, sigma33
	pk2Local[2] = hydroDiag + fiberPrefactor*Bmat[8];
	// in voigt notation, sigma23
	pk2Local[3] = fiberPrefactor*Bmat[7];
	// in voigt notation, sigma13
	pk2Local[4] = fiberPrefactor*Bmat[6];
	// in voigt notation, sigma12
	pk2Local[5] = fiberPrefactor*Bmat[3];
#ifdef DEBUG
  if(debug && 0){
    printf("Element ID = %d, gp = %d\n", e, gp);
    for(int i=0;i<6;i++){
      int index = pk2ptr[e]+6*gp+i;
      printf("PK2[%d] = %10.6e\n", index, pk2[index]);
    }
  }
#endif //DEBUG
  // TODO(Anil) : Move allocation and free to code initial setup
  free(Bmat);
}
