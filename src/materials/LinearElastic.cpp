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
    int pide = pid[e];
    double mu = properties[MAXMATPARAMS * pide + 1];
    double lambda = properties[MAXMATPARAMS * pide + 2];
    const unsigned int index = fptr[e] + ndim * ndim * gp;
    // Computation based on
    // http://run.usc.edu/femdefo/sifakis-courseNotes-TheoryAndDiscretization.pdf
    // section 3.2

    // Compute H matrix of the element from B^0 Matrix
    // Reuse mat1 for H computation
    double *H = mat1;
    ComputeH(e, gp, H);

    // Compute strain \epsilon = 0.5*(F+F^T)-I
    // This can be written as \epsilon = 0.5*(H+H^T)
    double *sigma_e = mat2; // Reuse mat2 for sigma computation
    // trace(\epsilon) = trace(F)-3 = trace(H)
    const double trEpsLambda = lambda*(H[0]+H[4]+H[8]);
    // \sigma = \lambda tr(\epsilon) I + 2 \mu \epsilon
    // \sigma = \lambda tr(H) I + \mu (H+H^T)
    // Compute 2*\mu*\epsilon
    for (int i = 0; i < ndim; ++i) {
      for (int j = 0; j < ndim; ++j) {
        const int indexL = j + i*ndim;
        sigma_e[indexL] = mu*(H[indexL]+H[i+j*ndim]);
      }
    }
    sigma_e[0] = sigma_e[0] + trEpsLambda;
    sigma_e[4] = sigma_e[4] + trEpsLambda;
    sigma_e[8] = sigma_e[8] + trEpsLambda;

    // Compute and store F = H + I
    double * const F_element_gp = &(F[index]);
    for (unsigned int i = 0; i < ndim2; ++i) {
      F_element_gp[i] = H[i];
    }
    F_element_gp[0] = F_element_gp[0] + 1.0;
    F_element_gp[4] = F_element_gp[4] + 1.0;
    F_element_gp[8] = F_element_gp[8] + 1.0;

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
