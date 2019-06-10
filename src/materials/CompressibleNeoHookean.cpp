#include "FemTech.h"
#include "blas.h"

// plane strain or three-dimensional compressible neo-Hookean
// Evaluates the Cauchy stress tensor

void CompressibleNeoHookean(int e, int gp){
//good example of how to reference F
	if (debug && 1==1) {
		printf("shp array e.%d with %d Gauss points, each with %d shp functions \n", e, GaussPoints[e], nShapeFunctions[e]);
		//printf("int.%d:\n", j);
		if(1==0){
			for (int k = 0; k < nShapeFunctions[e]; k++) {
				//printf("%8.5f ", shp[gptr[e] + j * GaussPoints[e] + k]);
				printf(" shp: %4.4f dshp: %8.4f %8.4f %8.4f\n",
						shp[gptr[e]   + gp * GaussPoints[e] + k],
						dshp[dsptr[e] + gp * GaussPoints[e] * ndim + k * ndim + 0],
						dshp[dsptr[e] + gp * GaussPoints[e] * ndim + k * ndim + 1],
						dshp[dsptr[e] + gp * GaussPoints[e] * ndim + k * ndim + 2]);
			}
		}
		//printf("\n");
		if(1==1){
			printf("Deformation Gradient, F for Gauss Point %d\n",gp);
			for(int i=0;i<ndim;i++){
					for(int j=0;j<ndim;j++){
						int index = fptr[e] + ndim*ndim*gp + ndim*i+j;
						printf(" F[%d]:%3.3e   ",index,F[index]);
					}
					printf("\n");
			}
			printf("\n");
		}
		int index2 = detFptr[e]+gp;
		printf("detF[%d]\n",index2);
	}

	//double J = 1;
	//printf("element %d, gauss point %d\n",e,gp);

	if(ndim == 2){
		// 6 values saved per gauss point for 3d
		for(int i=0;i<3;i++){
			int index = cptr[e]+3*gp+i;
			//cauchy[index] = 0.0
		  //printf("cauchy[%d]\n",index);
		}
	}
	if(ndim == 3){
    int index = fptr[e] + ndim * ndim * gp;
    int index2 = detFptr[e] + gp;
    double mu = properties[MAXMATPARAMS * e + 1];
    double lambda = properties[MAXMATPARAMS * e + 2];
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
		cauchy[cptr[e]+6*gp+0] = mu*(1.0-Cinv[0]) + logJ*Cinv[0];
			// in voigt notation, sigma22
		cauchy[cptr[e]+6*gp+1] = mu*(1.0-Cinv[4]) + logJ*Cinv[4];
			// in voigt notation, sigma33
		cauchy[cptr[e]+6*gp+2] = mu*(1.0-Cinv[8]) + logJ*Cinv[8];
			// in voigt notation, sigma23
		cauchy[cptr[e]+6*gp+3] = -mu*Cinv[7] + logJ*Cinv[7];
			// in voigt notation, sigma13
		cauchy[cptr[e]+6*gp+4] = -mu*Cinv[6] + logJ*Cinv[6];
			// in voigt notation, sigma12
		cauchy[cptr[e]+6*gp+5] = -mu*Cinv[3] + logJ*Cinv[3];

		for(int i=0;i<6;i++){
			int index = cptr[e]+6*gp+i;
			printf("cauchy[%d] = %3.3e\n",index,cauchy[index]);
		}
    free(C);
    free(Cinv);
	}
	return;
}
