#include "FemTech.h"
#include "blas.h"

void InternalForceUpdate(int e, int gp, double *force) {
  int cSize = 6;
  int nNodes = nShapeFunctions[e];
  double *fintGQ = (double*)calloc(nNodes*ndim, sizeof(double));
  int bColSize = nNodes*ndim;
  int Bsize = bColSize*cSize;
  double *B = (double*)calloc(Bsize, sizeof(double));
  // Calculate B matrix for each shape function
  for (int k = 0; k < nNodes; ++k) {
    StrainDisplacementMatrix(e, gp, k, &(B[6*ndim*k]));
  }
  double* sigma = &(cauchy[cptr[e]+6*gp]);
  // Compute B^T*sigma^n
  dgemv_(chy, &cSize, &bColSize, &one, B, &cSize, sigma, \
      &oneI, &zero, fintGQ, &oneI);
  // Add to fint
  int wIndex = gpPtr[e]+gp;
  // Following Belytschko equation 4.9.22
  const double preFactor = gaussWeights[wIndex]*detJacobian[wIndex];
  for (int k = 0; k < bColSize; ++k) {
    force[k] += preFactor*fintGQ[k];
  }
  printf("---- Internal Force Update ---\n");
  printf("Element : %d, GP : %d\n", e, gp);
  printf("B Matrix Transpose ---------\n");
  for (int j = 0; j < bColSize; ++j) {
    for (int i = 0; i < cSize; ++i) {
      printf("%12.6f  ", B[i+j*cSize]);
    }
    printf("\n");
  }
  for (int k = 0; k < bColSize; ++k) {
    printf("NodeNo : %d, dim = %d, Local force : %12.6f, Cum force : %12.6f\n", k/3, k%3, fintGQ[k], force[k]);
  }
  for (int k = 0; k < cSize; ++k) {
    printf("Start : %d, Stress value : %12.6f\n", cptr[e]+6*gp, sigma[k]);
  }
  printf("---- ----\n");
  free(fintGQ);
  free(B);
	return;
}
