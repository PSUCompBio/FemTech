#include "FemTech.h"
#include "blas.h"

// plane strain or three-dimensional Linear Elastic Material
// Evaluates the Cauchy stress tensor

void LinearElastic(int e, int gp) {
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
    int nNodes = nShapeFunctions[e];
    // double mu = properties[MAXMATPARAMS * e + 1];
    // double lambda = properties[MAXMATPARAMS * e + 2];
    int bColSize = nNodes*ndim;
    int cSize = 6;
    int Bsize = bColSize*cSize;
    // Variables to store intermediate outputs
    double *Bdn = (double*)calloc(cSize, sizeof(double));
    double *CBdn = (double*)calloc(cSize, sizeof(double));
    double *B = (double*)calloc(Bsize, sizeof(double));
    double *localDisplacement = (double*)calloc(nNodes*ndim, sizeof(double));
    for (int k = 0; k < nNodes; ++k) {
      int dIndex = connectivity[eptr[e]+k];
      localDisplacement[k*ndim+0] = displacements[dIndex*ndim+0];
      localDisplacement[k*ndim+1] = displacements[dIndex*ndim+1];
      localDisplacement[k*ndim+2] = displacements[dIndex*ndim+2];
    }
    // Calculate B matrix for each shape function
    for (int k = 0; k < nShapeFunctions[e]; ++k) {
      StressDisplacementMatrix(e, gp, k, &(B[6*ndim*k]));
    }
    // Compute B*d^n
    dgemv_(chn, &cSize, &bColSize, &one, B, &cSize, localDisplacement, \
        &oneI, &zero, Bdn, &oneI);
    // Compute sigma^n = C*B*d^n
    dgemv_(chn, &cSize, &cSize, &one, C, &cSize, Bdn, \
        &oneI, &zero, CBdn, &oneI);
		// 6 values saved per gauss point for 3d
		// in voigt notation, sigma11
		cauchy[cptr[e]+6*gp+0]=CBdn[0];
			// in voigt notation, sigma22
		cauchy[cptr[e]+6*gp+1]=CBdn[1];
			// in voigt notation, sigma33
		cauchy[cptr[e]+6*gp+2]=CBdn[2];
			// in voigt notation, sigma23
		cauchy[cptr[e]+6*gp+3]=CBdn[3];
			// in voigt notation, sigma13
		cauchy[cptr[e]+6*gp+4]=CBdn[4];
			// in voigt notation, sigma12
		cauchy[cptr[e]+6*gp+5]=CBdn[5];

    free(Bdn);
    free(CBdn);
    free(B);
    free(localDisplacement);

		// for(int i=0;i<6;i++){
		// 	int index = cptr[e]+6*gp+i;
		// 	printf("cauchy[%d] = %3.3e\n",index,cauchy[index]);
		// }
	}
	//printf("\n");
	return;
}
