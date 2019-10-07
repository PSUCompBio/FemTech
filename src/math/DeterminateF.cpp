#include "FemTech.h"
#include "blas.h"

void DeterminateF(int e, int gp){
#ifdef DEBUG
	if (debug && 1==0) {
			FILE_LOG_SINGLE(DEBUGLOG, "shp array e.%d with %d Gauss points, each with %d shp functions", e, GaussPoints[e], nShapeFunctions[e]);
			if(1==0){
				for (int k = 0; k < nShapeFunctions[e]; k++) {
					FILE_LOG_SINGLE(DEBUGLOG, "shp: %4.4f dshp: %8.4f %8.4f %8.4f",
							shp[gptr[e]   + gp * GaussPoints[e] + k],
							dshp[dsptr[e] + gp * GaussPoints[e] * ndim + k * ndim + 0],
							dshp[dsptr[e] + gp * GaussPoints[e] * ndim + k * ndim + 1],
							dshp[dsptr[e] + gp * GaussPoints[e] * ndim + k * ndim + 2]);
				}
			}
			if(1==0){
				FILE_LOG_SINGLE(DEBUGLOG, "Deformation Gradient, F for Gauss Point %d", gp);
				for(int i=0;i<ndim;i++){
						for(int j=0;j<ndim;j++){
							int index = fptr[e] + ndim*ndim*gp + ndim*i+j;
							printf(" F[%d]:%3.3e   ",index,F[index]);
						}
						printf("\n");
				}
				printf("\n");
			}
			int index2 = gpPtr[e]+gp;
			printf("detF[%d]\n",index2);
	}
#endif //DEBUG

  int index2 = gpPtr[e]+gp;

  if (ndim == 2) {
		int index = fptr[e] + ndim*ndim*gp;
		double da,db,dc,dd;
		//  [a b]
		//  [c d]
		//  det = (a*d-b*c)
		da=F[index + ndim*0+0]; //F11
		db=F[index + ndim*0+1]; //F12
		dc=F[index + ndim*1+0]; //F21
		dd=F[index + ndim*1+1]; //F22
    //detA = A[0],0]*A(1,1) - A(0,1)*A(1,0);
    //return detA;
		detF[index2] = da*dd-db*dc;
  }
  if (ndim == 3) {
		int index = fptr[e]+ndim*ndim*gp;
		double da,db,dc,dd,de,df,dg,dh,di;
		//  [a b c]
		//  [d e f]
		//  [g h i]
    //  |A| = a(ei − fh) − b(di − fg) + c(dh − eg)

		da=F[index + 0]; //F11
		db=F[index + 3]; //F12
		dc=F[index + 6]; //F13
		dd=F[index + 1]; //F21
		de=F[index + 4]; //F22
		df=F[index + 7]; //F23
		dg=F[index + 2]; //F31
		dh=F[index + 5]; //F32
		di=F[index + 8]; //F33

    detF[index2] = da*(de*di - df*dh) - db*(dd*di - df*dg) + dc*(dd*dh - de*dg);
  }
  else {
    FILE_LOG_SINGLE(ERROR, "ALERT: PROBLEM CALCULATING DETERMINANT OF MATRIX.\nMATRIX SIZE IS NOT 2X2 OR 3X3. SIMULATION CANCELLED");
    exit(0);
  }

	return;
}
