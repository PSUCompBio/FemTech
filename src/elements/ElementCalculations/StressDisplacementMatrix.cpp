#include "FemTech.h"

void StressDisplacementMatrix(int e,int gp){
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


	//FT = F.transpose();
	if(ndim == 2){
		// B in 2D is a 3 x 2 matrix
		double B11,B12,B21,B22,B31,B32;
		double dnIdx = 0.0;
		double dnIdy = 0.0;
		double F11 = 0.0;
		double F21 = 0.0;
		double F12 = 0.0;
		double F22 = 0.0;

		B11 = dnIdx * F11;
		B12 = dnIdx * F21;

		B21 = dnIdy * F12;
		B22 = dnIdy * F22;

		B31 = dnIdx * F12 + dnIdy * F11;
		B32 = dnIdx * F22 + dnIdy * F21;
	}
	if(ndim == 3){
		// B matrix in 3D is a 6 x 3 matrix
		double B11,B12,B13,B21,B22,B23,B31,B32,B33;
		double dnIdx = 0.0;
		double dnIdy = 0.0;
		double dnIdz = 0.0;

		double F11 = 0.0;
		double F21 = 0.0;
		double F31 = 0.0;

		double F12 = 0.0;
		double F22 = 0.0;
		double F32 = 0.0;

		double F13 = 0.0;
		double F23 = 0.0;
		double F33 = 0.0;

		int index = fptr[e] + ndim*ndim*gp;
		double F_local[3][3];

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

		for (int k = 0; k < nShapeFunctions[e]; k++) {
			//printf("%8.5f ", shp[gptr[e] + j * GaussPoints[e] + k]);
			printf(" shp: %4.4f dshp: %8.4f %8.4f %8.4f\n",
					shp[gptr[e]   + gp * GaussPoints[e] + k],
					dshp[dsptr[e] + gp * GaussPoints[e] * ndim + k * ndim + 0],
					dshp[dsptr[e] + gp * GaussPoints[e] * ndim + k * ndim + 1],
					dshp[dsptr[e] + gp * GaussPoints[e] * ndim + k * ndim + 2]);

			int index3=dsptr[e] + gp * GaussPoints[e] * ndim + k * ndim;
			double dnIdx = dshp[index3+0];
			double dnIdy = dshp[index3+1];
			double dnIdz = dshp[index3+2];
		}

		B11 = dnIdx * F11;
		B12 = dnIdx * F21;
		B13 = dnIdx * F31;

		B21 = dnIdy * F12;
		B22 = dnIdy * F22;
		B23 = dnIdy * F32;

		B31 = dnIdy*F13+dnIdz*F12 + dnIdy*F21+dnIdz*F22 + dnIdy*F33+dnIdz*F32;
		B32 = dnIdx*F13+dnIdz*F11 + dnIdx*F23+dnIdz*F21 + dnIdx*F33+dnIdz*F31;
		B33 = dnIdx*F12+dnIdy*F11 + dnIdx*F22+dnIdy*F21 + dnIdx*F32+dnIdy*F31;
	}

	return;
}
