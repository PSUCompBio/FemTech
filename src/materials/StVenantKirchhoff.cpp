#include "FemTech.h"
#include "blas.h"

// St. Venant-Kirchhoff
// Evaluates the Cauchy stress tensor

void StVenantKirchhoff(int e, int gp){

//good example of how to reference F
	if (debug && 1==0) {
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
    int index = fptr[e] + ndim*ndim*gp;
		int index2 = detFptr[e]+gp;
		double F_local[3][3];
    double FT_local[3][3];
		double b_local[3][3];
		double C_local[3][3];
		double E_local[3][3];
		double S_local[3][3];
		double FS_local[3][3];
		double FSFT_local[3][3];
		double cauchy_local[3][3];
		double mu = properties[MAXMATPARAMS*e+1]; // in future will be equal to component from  properties array
		double lambda = properties[MAXMATPARAMS*e+2]; // in future will be equal to component from  properties array
		double J=detF[index2];

		// local deformation gradient, F
    F_local[0][0]=F[index + ndim*0+0];
    F_local[0][1]=F[index + ndim*0+1];
		F_local[0][2]=F[index + ndim*0+2];
    F_local[1][0]=F[index + ndim*1+0];
		F_local[1][1]=F[index + ndim*1+1];
		F_local[1][2]=F[index + ndim*1+2];
		F_local[2][0]=F[index + ndim*2+0];
		F_local[2][1]=F[index + ndim*2+1];
		F_local[2][2]=F[index + ndim*2+2];

    // local transpose of deformation gradient, FT
		FT_local[0][0]=F[index + ndim*0+0];
 	 	FT_local[0][1]=F[index + ndim*1+0];
 	 	FT_local[0][2]=F[index + ndim*2+0];
 	 	FT_local[1][0]=F[index + ndim*0+1];
 	 	FT_local[1][1]=F[index + ndim*1+1];
 	 	FT_local[1][2]=F[index + ndim*2+1];
 	 	FT_local[2][0]=F[index + ndim*0+2];
 	 	FT_local[2][1]=F[index + ndim*1+2];
 	 	FT_local[2][2]=F[index + ndim*2+2];

		// C = F^T * F
		for (int i = 0; i < ndim; i++) {
			 for (int j = 0; j < ndim; j++) {
					double sum = 0.0;
					for (int k = 0; k < ndim; k++) {
						 sum = sum + FT_local[i][k] * F_local[k][j];
					}
					C_local[i][j] = sum;
			 }
		}

		//Compute Green-Lagrange Tensor: E= (1/2)*(F^T*F - I)
		E_local[0][0]=0.5*(C_local[0][0]-1.0);
		E_local[0][1]=0.5*(C_local[0][1]);
		E_local[0][2]=0.5*(C_local[0][2]);
		E_local[1][0]=0.5*(C_local[1][0]);
		E_local[1][1]=0.5*(C_local[1][1]-1.0);
		E_local[1][2]=0.5*(C_local[1][2]);
		E_local[2][0]=0.5*(C_local[2][0]);
		E_local[2][1]=0.5*(C_local[2][1]);
		E_local[2][2]=0.5*(C_local[2][2]-1.0);

		// Compute 2nd Piola-Kirchhoff Stress
		// lamda*tr(E) *I + 2*mu*E
		double trE=E_local[0][0]+E_local[1][1]+E_local[2][2];

		S_local[0][0]=lambda*trE + 2.0*mu*E_local[0][0];
		S_local[0][1]=             2.0*mu*E_local[0][1];
		S_local[0][2]=             2.0*mu*E_local[0][2];
		S_local[1][0]=             2.0*mu*E_local[1][0];
		S_local[1][1]=lambda*trE + 2.0*mu*E_local[1][1];
		S_local[1][2]=             2.0*mu*E_local[1][2];
		S_local[2][0]=             2.0*mu*E_local[2][0];
		S_local[2][1]=             2.0*mu*E_local[2][1];
		S_local[2][2]=lambda*trE + 2.0*mu*E_local[2][2];

    // sigma = 1/J * F * S * F^T
		// Compute F * S
	 for (int i = 0; i < ndim; i++) {
			 for (int j = 0; j < ndim; j++) {
					double sum = 0.0;
					for (int k = 0; k < ndim; k++) {
						 sum = sum + F_local[i][k] * S_local[k][j];
					}
					FS_local[i][j] = sum;
			 }
		}
		// Compute FS * F^T
	 for (int i = 0; i < ndim; i++) {
			 for (int j = 0; j < ndim; j++) {
					double sum = 0.0;
					for (int k = 0; k < ndim; k++) {
						 sum = sum + FS_local[i][k] * FT_local[k][j];
					}
					FSFT_local[i][j] = sum;
			 }
		}
		// 6 values saved per gauss point for 3d
		// in voigt notation, sigma11
		cauchy[cptr[e]+6*gp+0]=(1/J)*FSFT_local[0][0];
			// in voigt notation, sigma22
		cauchy[cptr[e]+6*gp+1]=(1/J)*FSFT_local[1][1];
			// in voigt notation, sigma33
		cauchy[cptr[e]+6*gp+2]=(1/J)*FSFT_local[2][2];
			// in voigt notation, sigma23
		cauchy[cptr[e]+6*gp+3]=(1/J)*FSFT_local[1][2];
			// in voigt notation, sigma13
		cauchy[cptr[e]+6*gp+4]=(1/J)*FSFT_local[0][2];
			// in voigt notation, sigma12
		cauchy[cptr[e]+6*gp+5]=(1/J)*FSFT_local[0][1];



		// for debugging can be removed , good example of how to reference F
		if (debug && 1==0){
			// print F
			for(int i=0;i<ndim;i++){
				for(int j=0;j<ndim;j++){
					printf(" F[%d%d]:%3.3e   ",i,j,F_local[i][j]);
				}
				printf("\n");
			}
			printf("\n");
			// Print F transpose
			for(int i=0;i<ndim;i++){
				for(int j=0;j<ndim;j++){
					printf(" FT[%d%d]:%3.3e   ",i,j,FT_local[i][j]);
				}
				printf("\n");
			}
			printf("\n");
			// Print right Cauchy Green tensor, C
			for(int i=0;i<ndim;i++){
				for(int j=0;j<ndim;j++){
					printf(" C[%d%d]:%3.3e   ",i,j,C_local[i][j]);
				}
				printf("\n");
			}
			printf("\n");

			// Print 2nd Piola-Kirchhoff stress tensor, S
			for(int i=0;i<ndim;i++){
				for(int j=0;j<ndim;j++){
					printf(" S[%d%d]:%3.3e   ",i,j,S_local[i][j]);
				}
				printf("\n");
			}
			printf("\n");

			// Print Cauchy stress tensor, C
			for(int i=0;i<ndim;i++){
				for(int j=0;j<ndim;j++){
					printf(" s[%d%d]:%3.3e   ",i,j,(1/J)*FSFT_local[i][j]);
				}
				printf("\n");
			}
			printf("\n");

		}

		if (debug && 1==0){
			for(int i=0;i<6;i++){
				int index = cptr[e]+6*gp+i;
				printf("cauchy[%d] = %3.3e\n",index,cauchy[index]);
			}
		}

	} // if ndim == 3

	return;


}
