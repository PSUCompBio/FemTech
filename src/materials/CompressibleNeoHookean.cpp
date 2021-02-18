#include "FemTech.h"
#include "blas.h"

// plane strain or three-dimensional compressible neo-Hookean
// Evaluates the PK2 stress tensor

void CompressibleNeoHookean(int e, int gp){
	if(ndim == 2){
		// 6 values saved per gauss point for 3d
		for(int i=0;i<3;i++){
			int index = pk2[e]+3*gp+i;
		}
	}
	if(ndim == 3){
    int pide = pid[e];
    double mu = properties[MAXMATPARAMS * pide + 1];
    double lambda = properties[MAXMATPARAMS * pide + 2];

		//From Bonet and Wood - Flagshyp
		//mu              = properties(2);
		//lambda          = properties(3);
		//J               = kinematics.J;
		//b               = kinematics.b;
		//Cauchy          = (mu/J)*(b - cons.I) + (lambda/J)*log(J)*cons.I;
    //C = F^T F
    // Bonet and Wood Equation 5.28
		//S          = mu*(I - C^{-1}) + lambda*ln(J)*C^{-1};

    // (b-I) = (H + H^T + H*H^T)
    double *H = mat1;
    ComputeH(e, gp, H);
    double *bmI = mat2;
    // Compute (H+H^T)
    for (int i = 0; i < ndim; ++i) {
      for (int j = 0; j < ndim; ++j) {
        const int indexL = j + i*ndim;
        bmI[indexL] = (H[indexL]+H[i+j*ndim]);
      }
    }
    // Compute H H^T and add to bmI
    dgemm_(chn, chy, &ndim, &ndim, &ndim, &one, H, &ndim,
           H, &ndim, &one, bmI, &ndim);
    // Store F in H
    // F = H + I
    H[0] = H[0] + 1.0;
    H[4] = H[4] + 1.0;
    H[8] = H[8] + 1.0;
    // J = det(F)
    const double J = det3x3Matrix(H);
    const double factor1 = mu/J;
    const double factor2 = lambda*log(J)/J;

    double *sigma_e = mat3;
    for (unsigned int i = 0; i < ndim2; ++i) {
      sigma_e[i] = factor1*bmI[i];
    }
    sigma_e[0] = sigma_e[0] + factor2;
    sigma_e[4] = sigma_e[4] + factor2;
    sigma_e[8] = sigma_e[8] + factor2;

		// 6 values saved per gauss point for 3d
		// in voigt notation, sigma11
    sigma_n[0] = sigma_e[0];
    // in voigt notation, sigma22
    sigma_n[1] = sigma_e[4];
    // in voigt notation, sigma33
    sigma_n[2] = sigma_e[8];
    // in voigt notation, sigma23
    sigma_n[3] = sigma_e[7];
    // in voigt notation, sigma13
    sigma_n[4] = sigma_e[6];
    // in voigt notation, sigma12
    sigma_n[5] = sigma_e[3];
	}
	return;
}
