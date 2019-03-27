#include "FemTech.h"
#include "blas.h"

// plane strain or three-dimensional compressible neo-Hookean
// Evaluates the Cauchy stress tensor

void CompressibleNeoHookean(int e, int gp){
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
		double cauchy_local[3][3];
		double mu = 100.0; // in future will be equal to component from  properties array
		double lambda = 166.6667; // in future will be equal to component from  properties array
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

    //MultiplyMatrices(FT_local,F_local,ndim,&d_local);
		for (int i = 0; i < ndim; i++) {
			 for (int j = 0; j < ndim; j++) {
					double sum = 0.0;
					for (int k = 0; k < ndim; k++) {
						 sum = sum + FT_local[i][k] * F_local[k][j];
					}
					b_local[i][j] = sum;
			 }
		}

		// for debugging can be removed , good example of how to reference F
		if (debug && 1==1){
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
			// Print left Cauchy Green tensor, b
			for(int i=0;i<ndim;i++){
				for(int j=0;j<ndim;j++){
					printf(" b[%d%d]:%3.4e   ",i,j,b_local[i][j]);
				}
				printf("\n");
			}
			printf("\n");
		}
		//printf("\n");

		//From Bonet and Wood - Flagshyp
		//mu              = properties(2);
		//lambda          = properties(3);
		//J               = kinematics.J;
		//b               = kinematics.b;
		//Cauchy          = (mu/J)*(b - cons.I) + (lambda/J)*log(J)*cons.I;

		// 6 values saved per gauss point for 3d
		// in voigt notation, sigma11
		cauchy[cptr[e]+6*gp+0]=(mu/J)*(b_local[0][0]-1.0) + (lambda/J)*log(J)*1.0;
			// in voigt notation, sigma22
		cauchy[cptr[e]+6*gp+1]=(mu/J)*(b_local[1][1]-1.0) + (lambda/J)*log(J)*1.0;
			// in voigt notation, sigma33
		cauchy[cptr[e]+6*gp+2]=(mu/J)*(b_local[2][2]-1.0) + (lambda/J)*log(J)*1.0;
			// in voigt notation, sigma23
		cauchy[cptr[e]+6*gp+3]=(mu/J)*(b_local[1][2]-0.0) + (lambda/J)*log(J)*0.0;
			// in voigt notation, sigma13
		cauchy[cptr[e]+6*gp+4]=(mu/J)*(b_local[0][2]-0.0) + (lambda/J)*log(J)*0.0;
			// in voigt notation, sigma12
		cauchy[cptr[e]+6*gp+5]=(mu/J)*(b_local[0][1]-0.0) + (lambda/J)*log(J)*0.0;

		for(int i=0;i<6;i++){
			int index = cptr[e]+6*gp+i;
			printf("cauchy[%d] = %3.4f\n",index,cauchy[index]);
		}
	}
	//printf("\n");
	return;


}
