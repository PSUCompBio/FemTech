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
    StressDisplacementMatrix(e, gp, k, &(B[6*ndim*k]));
  }
  double* sigma = &(cauchy[cptr[e]+6*gp]);
  // Compute B^T*sigma^n
  dgemv_(chy, &cSize, &bColSize, &one, B, &cSize, sigma, \
      &oneI, &zero, fintGQ, &oneI);
  // Add to fint
  int wIndex = gpPtr[e]+gp;
  const double preFactor = gaussWeights[wIndex]*detF[wIndex];
  for (int k = 0; k < bColSize; ++k) {
    force[k] += preFactor*fintGQ[k];
  }
  free(fintGQ);
	return;
}
