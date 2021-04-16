#include "FemTech.h"
#include "blas.h"

// plane strain or three-dimensional compressible neo-Hookean
// Evaluates the PK2 stress tensor

void CompressibleNeoHookeanAbaqus(int e, int gp){
	FILE *fp;
	fp=fopen("stress.txt","a");
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

		//From Abaqus
		//mu              = properties(2);
		//lambda          = properties(3);
		//J               = kinematics.J;
		//b               = kinematics.b;
    //B = F F^T
		//sigma          = (B - trace(B)*I/3)*mu/J^(5/3) + K*(J-1)*I;

    double *H = mat1;
    ComputeH(e, gp, H);
    double *bmI = mat2;
		double temp1, temp2, temp3;
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
		bmI[0] = bmI[0] + 1.0;
		bmI[4] = bmI[4] + 1.0;
		bmI[8] = bmI[8] + 1.0;
    // Store F in H
    // F = H + I
//		printf("%.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n", H[0], H[1], H[2], H[3], H[4], H[5], H[6], H[7], H[8], Time);

		H[0] = H[0] + 1.0;
    H[4] = H[4] + 1.0;
    H[8] = H[8] + 1.0;
	//	if(e==0)
	//		printf("%.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n", H[0], H[1], H[2], H[3], H[4], H[5], H[6], H[7], H[8], Time);

		temp1 = bmI[0] - (bmI[0]+bmI[4]+bmI[8])/3.0;
		temp2 = bmI[4] - (bmI[0]+bmI[4]+bmI[8])/3.0;
		temp3 = bmI[8] - (bmI[0]+bmI[4]+bmI[8])/3.0;
		bmI[0]=temp1;
		bmI[4]=temp2;
		bmI[8]=temp3;

	 // J = det(F)
    double J = det3x3Matrix(H);
    const double factor1 = mu/pow(J, 5.0/3.0);
    const double factor2 = (lambda+2*mu/3)*(J-1);

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
				fprintf(fp,"%.10f	%.10f %.10f %.10f %.10f %.10f %.10f\n", sigma_n[0], sigma_n[1], sigma_n[2], sigma_n[3], sigma_n[4], sigma_n[5], Time);
	//fprintf(fp, "%.10f %.10f\n", J , dt);
				}
						fclose(fp);
	return;
}
