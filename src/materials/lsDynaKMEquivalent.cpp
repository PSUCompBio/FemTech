#include "FemTech.h"
#include "blas.h"

/* three-dimensional compressible neo-Hookean model with one term prony series
 * equivalent to LS-Dyna implementation of KM model (MAT_061) with FO=0, ie Maxwell model.
 * G(t) = G_\infity + (G_0 - G_\infity) exp (-\beta t)
 * In prony series equivalent G_{eq} *(1 + g_1 exp (-\beta t))
 * g_1 = (G_0-G_\infity)/G_\infity, G_{eq} = G_\infity
 * Evaluates the cauchy stress tensor
 * Viscoelasticity implementation based on : 
 * Kaliske, M. and Rothert, H., 1997. Formulation and 
 * implementation of three-dimensional viscoelasticity at small and finite 
 * strains. Computational Mechanics, 19(3), pp.228-239. 
 * Expected material properties :
 * \rho      = properties(0)
 * \lambda   = properties(1)
 * G_\infity = properties(2)
 * g_1       = properties(3)
 * \beta_1   = properties(4)
 * */

void lsDynaKMEquivalent(int e, int gp) {
	if(ndim == 2) {
		// 6 values saved per gauss point for 3d
		for(int i=0;i<3;i++){
			int index = pk2[e]+3*gp+i;
		}
	}
	if(ndim == 3) {
    const unsigned int pideOffset = MAXMATPARAMS*pid[e];
    const double mu = properties[pideOffset + 2];
    const double lambda = properties[pideOffset + 1];
    const double g_1 = properties[pideOffset + 3];
    const double beta_1 = properties[pideOffset + 4];

    const unsigned int index = fptr[e] + ndim * ndim * gp;

		//From Bonet and Wood - Flagshyp
		//Cauchy          = (mu/J)*(b - cons.I) + (lambda/J)*log(J)*cons.I;
    // Bonet and Wood Equation 5.29
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
    // Compute and store F = H + I
    double * const F_element_gp = &(F[index]);
    for (unsigned int i = 0; i < ndim2; ++i) {
      F_element_gp[i] = H[i];
    }
    F_element_gp[0] = F_element_gp[0] + 1.0;
    F_element_gp[4] = F_element_gp[4] + 1.0;
    F_element_gp[8] = F_element_gp[8] + 1.0;
    // J = det(F)
    const double J = det3x3Matrix(F_element_gp);
    const double factor1 = mu/J;
    const double factor2 = lambda*log(J)/J;

    double *sigma_e = mat3;
    for (unsigned int i = 0; i < ndim2; ++i) {
      sigma_e[i] = factor1*bmI[i];
    }
    sigma_e[0] = sigma_e[0] + factor2;
    sigma_e[4] = sigma_e[4] + factor2;
    sigma_e[8] = sigma_e[8] + factor2;

    // minus of pressure computation
    const double mPressure = (sigma_e[0] + sigma_e[4] + sigma_e[8])/3.0;

    // Reuse temp storage mat1 to store Deviatoric part of sigma
    double *sigma_eDev = mat1;
    for (unsigned int i = 0; i < ndim2; ++i) {
      sigma_eDev[i] = sigma_e[i];
    }
    sigma_eDev[0] = sigma_eDev[0] - mPressure;
    sigma_eDev[4] = sigma_eDev[4] - mPressure;
    sigma_eDev[8] = sigma_eDev[8] - mPressure;

    // Access histroy dependence
    double *elemS = S0n[e];
    double *S0nLocal = &(elemS[gp*ndim2]);
    double *elemH = Hn[e];
    double *HnI = &(elemH[gp*ndim2]);
    // Compute C1i = exp(-dt/t_i) and C2i = g_i (1-exp(-dt/t_i))/(dt/t_i)
    const double rt = dt*beta_1;
    const double c1i = exp(-rt);
    const double c2i = g_1*(1.0-c1i)/rt;
    // Update H_j^{n+1} = c1j*H_j^n + c2j*[Dev S_0^{n+1} - S_0^n]
    // Update S^{n+1} = S_0^{n+1} + \sum_j H_j^{n+1}
    for (int j = 0; j < ndim2; ++j) {
      HnI[j] = c1i*HnI[j]+c2i*(sigma_eDev[j]-S0nLocal[j]);
      sigma_e[j] = sigma_e[j] + HnI[j];
    }

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

    // Update deviatoric part of S_0 for next time step
    for (int i = 0; i < ndim2; ++i) {
      S0nLocal[i] = sigma_eDev[i];
    }
	}
	return;
}
