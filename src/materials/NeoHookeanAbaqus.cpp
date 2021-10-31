#include "FemTech.h"
#include "blas.h"

// plane strain or three-dimensional compressible neo-Hookean
// Evaluates the Cauchy stress tensor

void NeoHookeanAbaqus(int e, int gp){
	if(ndim == 3) {
    int index = fptr[e] + ndim * ndim * gp;
    int index2 = detFptr[e] + gp;
    int pide = pid[e];
    double mu = properties[MAXMATPARAMS * pide + 1];
    double lambda = properties[MAXMATPARAMS * pide + 2];
    double J = detF[index2];

    // Compute Green-Lagrange Tensor: C = F^T*F
    double *Cmat = mat1;
    double *Cinv = mat2;
    double *F_element_gp = &(F[index]);
    dgemm_(chy, chn, &ndim, &ndim, &ndim, &one, F_element_gp, &ndim,
           F_element_gp, &ndim, &zero, Cmat, &ndim);
    double Cdet;
    inverse3x3Matrix(Cmat, Cinv, &Cdet);

		//S = \mu*J^{-2/3}*I + (-\mu*J^{-2/3}*tr(C)/3 + 
    // (\lambda+2*\mu/3)*(J-1))J*C^{-1};
    const double Jm23 = pow(J, -2.0/3.0);
    const double traceC = Cmat[0] + Cmat[4] + Cmat[8];
    const double preFactor = -mu*Jm23*traceC/3.0 + (lambda+2.0*mu/3.0)*(J-1.0)*J;

		// 6 values saved per gauss point for 3d
    // Computing PK2
		// in voigt notation, sigma11
		pk2[pk2ptr[e]+6*gp+0] = mu*Jm23+preFactor*Cinv[0];
			// in voigt notation, sigma22
		pk2[pk2ptr[e]+6*gp+1] = mu*Jm23+preFactor*Cinv[4];
			// in voigt notation, sigma33
		pk2[pk2ptr[e]+6*gp+2] = mu*Jm23+preFactor*Cinv[8];
			// in voigt notation, sigma23
		pk2[pk2ptr[e]+6*gp+3] = preFactor*Cinv[7];
			// in voigt notation, sigma13
		pk2[pk2ptr[e]+6*gp+4] = preFactor*Cinv[6];
			// in voigt notation, sigma12
		pk2[pk2ptr[e]+6*gp+5] = preFactor*Cinv[3];
    FILE_LOGArraySingle(DEBUGLOGIGNORE, &(pk2[pk2ptr[e]+6*gp]), 6, \
        "Element : %d, GP : %d Stress values", e, gp);
	}
	return;
}
