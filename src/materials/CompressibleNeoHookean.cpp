#include "FemTech.h"
#include "blas.h"

// plane strain or three-dimensional compressible neo-Hookean
// Evaluates the Cauchy stress tensor

void CompressibleNeoHookean(int e, int gp){
	if(ndim == 3){
    int pide = pid[e];
    double mu = properties[MAXMATPARAMS * pide + 1];
    double lambda = properties[MAXMATPARAMS * pide + 2];
    const unsigned int index = fptr[e] + ndim * ndim * gp;

		//From Bonet and Wood - Flagshyp
		//mu              = properties(2);
		//lambda          = properties(3);
		//Cauchy          = (mu/J)*(b - cons.I) + (lambda/J)*log(J)*cons.I;
    //C = F^T F
    // Bonet and Wood Equation 5.28
		//S          = mu*(I - C^{-1}) + lambda*ln(J)*C^{-1};

    // (b-I) = (H + H^T + H*H^T)
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

    double *b = mat2;
    // Compute F F^T and store to b 
    dgemm_(chn, chy, &ndim, &ndim, &ndim, &one, F_element_gp, &ndim,
           F_element_gp, &ndim, &zero, b, &ndim);
    // Compute b - I
    b[0] = b[0] - 1.0;
    b[4] = b[4] - 1.0;
    b[8] = b[8] - 1.0;

    // J = det(F)
    const double J = det3x3Matrix(F_element_gp);
    const double factor1 = mu/J;
    const double factor2 = lambda*log(J)/J;

    double *sigma_e = mat3;
    for (unsigned int i = 0; i < ndim2; ++i) {
      sigma_e[i] = factor1*b[i];
    }
    sigma_e[0] = sigma_e[0] + factor2;
    sigma_e[4] = sigma_e[4] + factor2;
    sigma_e[8] = sigma_e[8] + factor2;

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
	}
	return;
}
