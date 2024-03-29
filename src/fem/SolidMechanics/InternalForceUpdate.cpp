#include "FemTech.h"
#include "blas.h"
// Drupal : Check if the operating system is Windows
#ifdef _WIN32
#include "matrix_operations.h"
#endif

void InternalForceUpdate(int e, int gp, double *force) {
  int cSize = 6;
  int nNodesL = nShapeFunctions[e];
  int bColSize = nNodesL*ndim;
  int Bsize = bColSize*cSize;
  for (int k = 0; k < Bsize; ++k) {
    B[k] = 0.0;
  }
  // Calculate B matrix for each shape function
  for (int k = 0; k < nNodesL; ++k) {
    StrainDisplacementMatrix(e, gp, k, &(B[6*ndim*k]));
  }
  double* sigma = &(sigma_n[sigmaptr[e]+6*gp]);
  // Compute B^T*sigma^n
// Drupal : Check if the operating system is Windows
#ifdef _WIN32
  FEMTECH_dgemv_(chy, cSize, bColSize, one, B, cSize, sigma, oneI, zero, fintGQ, oneI);
#else
  dgemv_(chy, &cSize, &bColSize, &one, B, &cSize, sigma, &oneI, &zero, fintGQ, &oneI);
#endif

  
  // Add to fint
  int wIndex = gpPtr[e]+gp;
  // Following Belytschko equation 4.9.22
  const double preFactor = gaussWeights[wIndex]*detJacobian[wIndex];
  for (int k = 0; k < bColSize; ++k) {
    force[k] += preFactor*fintGQ[k];
  }
	return;
}
