#include "FemTech.h"
#include "blas.h"

void CalculateDeformationGradient(int e, int gp){

	//good example of how to reference deformation gradient, F
	if (debug && 1==0) {
			printf("shp array e.%d with %d Gauss points, each with %d shp functions \n", e, GaussPoints[e], nShapeFunctions[e]);
			//printf("int.%d:\n", j);
			for (int k = 0; k < nShapeFunctions[e]; k++) {
				//printf("%8.5f ", shp[gptr[e] + j * GaussPoints[e] + k]);
				printf(" shp: %4.4f dshp: %8.4f %8.4f %8.4f\n",
						shp[gptr[e]   + gp * GaussPoints[e] + k],
						dshp[dsptr[e] + gp * GaussPoints[e] * ndim + k * ndim + 0],
						dshp[dsptr[e] + gp * GaussPoints[e] * ndim + k * ndim + 1],
						dshp[dsptr[e] + gp * GaussPoints[e] * ndim + k * ndim + 2]);
			}
			//printf("\n");
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

	// following Bonet and Wood; F = xlocal*DN_X in flagshyp

    //printf("size of Xlocal columns: %d\n",nShapeFunctions[e]);
    //double xlocal[ndim][nShapeFunctions[e]];
    //double ds[nShapeFunctions[e]][ndim];
    //double defgrd[ndim][ndim];
    if( (eptr[e+1]-eptr[e]) != nShapeFunctions[e]){
				printf("Warning from CalculateDeformationGradient.cpp: (eptr[e+1]-eptr[e]) != nShapeFunctions[e])" );
        printf("Check it out, probally a bug in allocation" );
		}
		//int xlcounter = 0;
		//for (int a = eptr[e]; a < eptr[e + 1]; a++) {
		//		int node_a = connectivity[a];
			  //printf("d %d->%3.3e %3.3e %3.3e\n",a,coordinates[ndim*node_a+0]+displacements[ndim*node_a+0],
				//													   	  coordinates[ndim*node_a+1]+displacements[ndim*node_a+1],
				//													      coordinates[ndim*node_a+2]+displacements[ndim*node_a+2]);
		//		xlocal[0][xlcounter] = coordinates[ndim*node_a+0]+displacements[ndim*node_a+0];
		//		xlocal[1][xlcounter] = coordinates[ndim*node_a+1]+displacements[ndim*node_a+1];
		//		xlocal[2][xlcounter] = coordinates[ndim*node_a+2]+displacements[ndim*node_a+2];
				//printf("%3.3e %3.3e %3.3e\n", Xlocal[0][xlcounter],Xlocal[1][xlcounter],Xlocal[2][xlcounter] );
		//		xlcounter++;
   // } // loop on a


		//for (int i=0;i<nShapeFunctions[e];i++){
		//		for(int j=0;j<ndim;j++){
		//			int index2= dsptr[e] + gp * GaussPoints[e] * ndim + i * ndim + j;
    //      ds[i][j]=dshp[index2];
		//	}
	//	}
//

	//	for(int i=0;i<ndim;i++){
	//		for(int j=0;j<ndim;j++){
	//				defgrd[i][j]=0.0;
	//		}
	//	}

    double theSum = 0.0;
    for(int i=0;i<ndim;i++){
			for(int j=0;j<ndim;j++){
				theSum = 0.0;
				int index = fptr[e] + ndim*ndim*gp + ndim*i+j;
				for(int k=0;k<nShapeFunctions[e];k++){
						int a = eptr[e];
						int node_a = connectivity[a+k];
					//	theSum = theSum + xlocal[i][k] * ds[k][j];
						int index2= dsptr[e] + gp * GaussPoints[e] * ndim + k * ndim + j;
						//if(theSum < 1e-16) {theSum = 0.0;}
						//printf("1) %3.3e %3.3e\n",xlocal[i][k],coordinates[ndim*node_a+i]+displacements[ndim*node_a+i]);
						//printf("2) %3.3e %3.3e\n",ds[k][j],dshp[index2]);
						//int index = fptr[e] + ndim*ndim*gp + ndim*i+j;
						theSum = theSum +  (coordinates[ndim*node_a+i]+displacements[ndim*node_a+i]) * dshp[index2];
				} // loop on k
		//		defgrd[i][j]=theSum;
				F[index] = theSum;

			}
		}

	//	for(int i=0;i<ndim;i++){
	//		printf("%3.3e %3.3e %3.3e\n", defgrd[i][0],defgrd[i][1],defgrd[i][2]);
		//}

		//for(int i=0;i<ndim;i++){
		//	int index = fptr[e] + ndim*ndim*gp + ndim*i;
		//	printf(">%3.3e %3.3e %3.3e\n", F[index+0],F[index+1],F[index+2]);
		//}

  //printf("\n");

	if(debug && 1==0){
		printf("--------F for gauss point %d --------\n",gp);
		for(int i=0;i<ndim;i++){
			for(int j=0;j<ndim;j++){
					int index = fptr[e] + ndim*ndim*gp + ndim*i+j;
					printf("%3.3e   ",F[index]);
			} //loop on j
			printf("\n");
		} //loop on i
	} // if debug


	//free(ds);
  //free(Xlocal);
	return;
}
