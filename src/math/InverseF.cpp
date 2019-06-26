#include "FemTech.h"
#include "blas.h"

void InverseF(int e, int gp){
	//good example of how to reference F
#ifdef DEBUG
	if (debug && 1==0) {
			printf("shp array e.%d with %d Gauss points, each with %d shp functions \n", e, GaussPoints[e], nShapeFunctions[e]);
			//printf("int.%d:\n", j);
			if(1==1){
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
			if(1==0){
				printf("Inverse of Deformation Gradient, F for Gauss Point %d\n",gp);
				for(int i=0;i<ndim;i++){
						for(int j=0;j<ndim;j++){
							int index = fptr[e] + ndim*ndim*gp + ndim*i+j;
							printf(" invF[%d]:%3.3e   ",index,invF[index]);
						}
						printf("\n");
				}
				printf("\n");
			}
	}
#endif //DEBUG

	if(ndim == 2) {
		int index = fptr[e] + ndim*ndim*gp;
		int index2 = detFptr[e]+gp;
		double da,db,dc,dd;
		//  [a b]
		//  [c d]
    //  inv A = 1/detA * [d  -b]
		//  						   [-c  a]
		da = F[index + ndim*0+0];
		db = F[index + ndim*0+1];
		dc = F[index + ndim*1+0];
		dd = F[index + ndim*1+1];

		invF[index + ndim*0+0] = 1.0/detF[index2]*dd;
		invF[index + ndim*0+1] = -1.0*1.0/detF[index2]*db;
		invF[index + ndim*1+0] = 1.0*1.0/detF[index2]*dc;
		invF[index + ndim*1+1] = 1.0/detF[index2]*da;


	}
	if(ndim == 3){
		int index = fptr[e] + ndim*ndim*gp;
		//double da,db,dc,dd,de,df,dg,dh,di;
		double A[3][3];
		/*da=F[index + ndim*0+0]; //F11
		db=F[index + ndim*0+1]; //F12
		dc=F[index + ndim*0+2]; //F13
		dd=F[index + ndim*1+0]; //F21
		de=F[index + ndim*1+1]; //F22
		df=F[index + ndim*1+2]; //F23
		dg=F[index + ndim*2+0]; //F31
		dh=F[index + ndim*2+1]; //F32
		di=F[index + ndim*3+2]; //F33 */
		A[0][0]=F[index + ndim*0+0]; //F11
		A[0][1]=F[index + ndim*0+1]; //F12
		A[0][2]=F[index + ndim*0+2]; //F13

		A[1][0]=F[index + ndim*1+0]; //F21
		A[1][1]=F[index + ndim*1+1]; //F22
		A[1][2]=F[index + ndim*1+2]; //F23

		A[2][0]=F[index + ndim*2+0]; //F31
		A[2][1]=F[index + ndim*2+1]; //F32
		A[2][2]=F[index + ndim*2+2]; //F33

		//     0  1  2
		// 0  [a  b  c]
		// 1  [d  e  f]
		// 2  [g  h  i]
    double detA = detF[detFptr[e]+gp];

		invF[index + ndim*0+0] = 1.0/detA*(A[1][1]*A[2][2] - A[1][2]*A[2][1]);
		invF[index + ndim*0+1] = 1.0/detA*(A[0][2]*A[2][1] - A[0][1]*A[2][2]);
		invF[index + ndim*0+2] = 1.0/detA*(A[0][1]*A[1][2] - A[0][2]*A[1][1]);

		invF[index + ndim*1+0] = 1.0/detA*(A[1][2]*A[2][0] - A[1][0]*A[2][2]);
		invF[index + ndim*1+1] = 1.0/detA*(A[0][0]*A[2][2] - A[0][2]*A[2][0]);
		invF[index + ndim*1+2] = 1.0/detA*(A[0][2]*A[1][0] - A[0][0]*A[1][2]);

		invF[index + ndim*2+0] = 1.0/detA*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);
		invF[index + ndim*2+1]=  1.0/detA*(A[0][1]*A[2][0] - A[0][0]*A[2][1]);
		invF[index + ndim*2+2] = 1.0/detA*(A[0][0]*A[1][1] - A[0][1]*A[1][0]);
	}

	return;
}
