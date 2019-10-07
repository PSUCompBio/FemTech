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
    int index = fptr[e] + ndim * ndim * gp;
    int index2 = detFptr[e] + gp;
    int pide = pid[e];
    double mu = properties[MAXMATPARAMS * pide + 1];
    double lambda = properties[MAXMATPARAMS * pide + 2];
    double J = detF[index2];

    // Compute Green-Lagrange Tensor: C = F^T*F
    double matSize = ndim * ndim;
    double *C = (double *)malloc(matSize * sizeof(double));
    double *Cinv = (double *)malloc(matSize * sizeof(double));
    double *F_element_gp = &(F[index]);
    dgemm_(chy, chn, &ndim, &ndim, &ndim, &one, F_element_gp, &ndim,
           F_element_gp, &ndim, &zero, C, &ndim);
    double Cdet;
    inverse3x3Matrix(C, Cinv, &Cdet);

		//From Bonet and Wood - Flagshyp
		//mu              = properties(2);
		//lambda          = properties(3);
		//J               = kinematics.J;
		//b               = kinematics.b;
		//Cauchy          = (mu/J)*(b - cons.I) + (lambda/J)*log(J)*cons.I;
    //C = F^T F
    // Bonet and Wood Equation 5.28
		//S          = mu*(I - C^{-1}) + lambda*ln(J)*C^{-1};

		// 6 values saved per gauss point for 3d
    // Computing PK2
    double logJ = lambda*log(J);
		// in voigt notation, sigma11
		pk2[pk2ptr[e]+6*gp+0] = mu*(1.0-Cinv[0]) + logJ*Cinv[0];
			// in voigt notation, sigma22
		pk2[pk2ptr[e]+6*gp+1] = mu*(1.0-Cinv[4]) + logJ*Cinv[4];
			// in voigt notation, sigma33
		pk2[pk2ptr[e]+6*gp+2] = mu*(1.0-Cinv[8]) + logJ*Cinv[8];
			// in voigt notation, sigma23
		pk2[pk2ptr[e]+6*gp+3] = -mu*Cinv[7] + logJ*Cinv[7];
			// in voigt notation, sigma13
		pk2[pk2ptr[e]+6*gp+4] = -mu*Cinv[6] + logJ*Cinv[6];
			// in voigt notation, sigma12
		pk2[pk2ptr[e]+6*gp+5] = -mu*Cinv[3] + logJ*Cinv[3];
    FILE_LOGArraySingle(DEBUGLOGIGNORE, &(pk2[pk2ptr[e]+6*gp]), 6, \
        "Element : %d, GP : %d Stress values", e, gp);
    free(C);
    free(Cinv);
	}
	return;
}
