#include "FemTech.h"
#include "blas.h"

// plane strain or three-dimensional Linear Elastic Material
// Evaluates the PK2 stress tensor

void LinearElastic(int e, int gp) {
	if(ndim == 2){
		// 6 values saved per gauss point for 3d
		for(int i=0;i<3;i++){
			int index = pk2ptr[e]+3*gp+i;
		}
	}
	if(ndim == 3){
    int index = fptr[e] + ndim * ndim * gp;
    int index2 = detFptr[e] + gp;
    int pide = pid[e];
    double mu = properties[MAXMATPARAMS * pide + 1];    
    double lambda = properties[MAXMATPARAMS * pide + 2];
    double J = detF[index2];

    // Compute strain \epsilon = 0.5*(F+F^T)-I
    double matSize = ndim * ndim;
    double *eps = (double *)malloc(matSize * sizeof(double));
    double *F_element_gp = &(F[index]);
    const double trEps = F_element_gp[0]+F_element_gp[4]+F_element_gp[8]-3.0;
    for (int i = 0; i < ndim; ++i) {
      for (int j = 0; j < ndim; ++j) {
        const int index = j + i*ndim;
        eps[index] = 0.5*(F_element_gp[index]+F_element_gp[i+j*ndim]);
      }
    }
    eps[0] -= 1.0;
    eps[4] -= 1.0;
    eps[8] -= 1.0;

    // Compute sigma = \lambda tr(\eps) I + 2 \mu \eps
    double *sigma = (double *)malloc(matSize * sizeof(double));
    for (int i = 0; i < matSize; ++i) {
      sigma[i] = 2.0 * mu * eps[i];
    }
    sigma[0] += lambda * trEps;
    sigma[4] += lambda * trEps;
    sigma[8] += lambda * trEps;

    // Compute pk2 : S = det(F) F^{-1} \sigma F^{-T}
    double *fInv = (double *)malloc(matSize * sizeof(double));
    double *S1 = (double *)malloc(matSize * sizeof(double));
    double *S = (double *)malloc(matSize * sizeof(double));
    InverseF(e, gp, fInv);
    // Compute F^{-1}*\sigma
    dgemm_(chn, chn, &ndim, &ndim, &ndim, &one, fInv, &ndim,
           sigma, &ndim, &zero, S1, &ndim);
    // Compute det(F)*F^{-1}*\sigma*F^{-T}
    dgemm_(chn, chy, &ndim, &ndim, &ndim, &J, S1, &ndim,
           fInv, &ndim, &zero, S, &ndim);
		// 6 values saved per gauss point for 3d
		// in voigt notation, sigma11
    pk2[pk2ptr[e] + 6 * gp + 0] = S[0];
    // in voigt notation, sigma22
    pk2[pk2ptr[e] + 6 * gp + 1] = S[4];
    // in voigt notation, sigma33
    pk2[pk2ptr[e] + 6 * gp + 2] = S[8];
    // in voigt notation, sigma23
    pk2[pk2ptr[e] + 6 * gp + 3] = S[7];
    // in voigt notation, sigma13
    pk2[pk2ptr[e] + 6 * gp + 4] = S[6];
    // in voigt notation, sigma12
    pk2[pk2ptr[e] + 6 * gp + 5] = S[3];

    free(eps);
    free(sigma);
    free(fInv);
    free(S);
	}
	return;
}
