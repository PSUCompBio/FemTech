#include "FemTech.h"
#include "blas.h"

//<<<<<<< addtrussV1
//void CalculateDeformationGradient(int e, int gp){

	//good example of how to reference deformation gradient, F
	//if (debug && 1==0) {
	//		printf("shp array e.%d with %d Gauss points, each with %d shp functions \n", e, GaussPoints[e], nShapeFunctions[e]);
			//printf("int.%d:\n", j);
	//		for (int k = 0; k < nShapeFunctions[e]; k++) {
				//printf("%8.5f ", shp[gptr[e] + j * GaussPoints[e] + k]);
	//			printf(" shp: %4.4f dshp: %8.4f %8.4f %8.4f\n",
	//					shp[gptr[e]   + gp * GaussPoints[e] + k],
	//					dshp[dsptr[e] + gp * GaussPoints[e] * ndim + k * ndim + 0],
	//					dshp[dsptr[e] + gp * GaussPoints[e] * ndim + k * ndim + 1],
	//					dshp[dsptr[e] + gp * GaussPoints[e] * ndim + k * ndim + 2]);
	//		}
	//		//printf("\n");
	//		printf("Deformation Gradient, F for Gauss Point %d\n",gp);
	//		for(int i=0;i<ndim;i++){
	//				for(int j=0;j<ndim;j++){
		//				int index = fptr[e] + ndim*ndim*gp + ndim*i+j;
		//				printf(" F[%d]:%3.3e   ",index,F[index]);
	//				}
	//				printf("\n");
	//		}
	//		printf("\n");
	//}

	// following Bonet and Wood; F = xlocal*DN_X in flagshyp

    //printf("size of Xlocal columns: %d\n",nShapeFunctions[e]);
    //double xlocal[ndim][nShapeFunctions[e]];
    //double ds[nShapeFunctions[e]][ndim];
    //double defgrd[ndim][ndim];
 // if( (eptr[e+1]-eptr[e]) != nShapeFunctions[e]){
	///		printf("Warning from CalculateDeformationGradient.cpp: (eptr[e+1]-eptr[e]) != nShapeFunctions[e])" );
    //  printf("Check it out, probally a bug in allocation" );
//	}

//  double theSum = 0.0;
 // for(int i=0;i<ndim;i++){
	//	for(int j=0;j<ndim;j++){
	//		theSum = 0.0;
//			int index = fptr[e] + ndim*ndim*gp + ndim*j+i;
//			for(int k=0;k<nShapeFunctions[e];k++){
//					int a = eptr[e];
//					int node_a = connectivity[a+k];
//					int index2= dsptr[e] + gp * GaussPoints[e] * ndim + k * ndim + j;
//					theSum = theSum +  (coordinates[ndim*node_a+i]+displacements[ndim*n/ode_a+i]) * dshp[index2];
//			} // loop on k
//			F[index] = theSum;

//		}
//	}




//	if(debug && 1==0){
//		printf("--------F for gauss point %d --------\n",gp);
//		for(int i=0;i<ndim;i++){
//			for(int j=0;j<ndim;j++){
//					int index = fptr[e] + ndim*ndim*gp + ndim*j+i;
//					printf("%3.3e   ",F[index]);
//			} //loop on j
//			printf("\n");
//		} //loop on i
//	} // if debug

//	return;
//=======
void CalculateDeformationGradient(int e, int gp) {

  // good example of how to reference deformation gradient, F
  if (debug && 1 == 0) {
    printf("shp array e.%d with %d Gauss points, each with %d shp functions \n",
           e, GaussPoints[e], nShapeFunctions[e]);
    for (int k = 0; k < nShapeFunctions[e]; k++) {
      printf(" shp: %4.4f dshp: %8.4f %8.4f %8.4f\n",
             shp[gptr[e] + gp * GaussPoints[e] + k],
             dshp[dsptr[e] + gp * GaussPoints[e] * ndim + k * ndim + 0],
             dshp[dsptr[e] + gp * GaussPoints[e] * ndim + k * ndim + 1],
             dshp[dsptr[e] + gp * GaussPoints[e] * ndim + k * ndim + 2]);
    }
  }

  // following Bonet and Wood; F = xlocal*DN_X in flagshyp
  if ((eptr[e + 1] - eptr[e]) != nShapeFunctions[e]) {
    printf("Warning from CalculateDeformationGradient.cpp: (eptr[e+1]-eptr[e]) "
           "!= nShapeFunctions[e])");
    printf("Check it out, probally a bug in allocation");
  }
  double theSum = 0.0;
  for (int i = 0; i < ndim; i++) {
    for (int j = 0; j < ndim; j++) {
      theSum = 0.0;
      int index = fptr[e] + ndim * ndim * gp + ndim * j + i;
      for (int k = 0; k < nShapeFunctions[e]; k++) {
        int a = eptr[e];
        int node_a = connectivity[a + k];
        int index2 = dsptr[e] + gp * GaussPoints[e] * ndim + k * ndim + i;
        theSum = theSum + (coordinates[ndim * node_a + j] +
                           displacements[ndim * node_a + j]) *
                              dshp[index2];
      } // loop on k
      F[index] = theSum;
    }
  }

  if (debug && 1 == 0) {
    printf("--------F for gauss point %d --------\n", gp);
    for (int i = 0; i < ndim; i++) {
      for (int j = 0; j < ndim; j++) {
        int index = fptr[e] + ndim * ndim * gp + ndim * j + i;
        printf("%3.3e   ", F[index]);
      } // loop on j
      printf("\n");
    } // loop on i
  }   // if debug

  return;
//>>>>>>> master
}
