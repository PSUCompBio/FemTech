#include "FemTech.h"

void InternalForceUpdate(int e, int gp){

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
	//StressDisplacementMatrix(i,j);
//f_int_e = f_int_e + ((disp_mat.transpose()) * sigma_e * wtx_normal[ix] * wty_normal[ix][iy] * wtz_normal[ix][iy][iz] * det_store[i][ix][iy][iz]);
	return;
}
