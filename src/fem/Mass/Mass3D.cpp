#include "FemTech.h"

#include "blas.h"

void MassElementMatrix(double *Me, int e) {
  // number of shape functions * ndim
  int bColSize = nShapeFunctions[e]*ndim;
  int mLocalSize = bColSize*bColSize;

  //TODO(Anil) remove allocations inside the loop
  double *MeGQ = (double*)calloc(mLocalSize, sizeof(double));
  int nColSize = nShapeFunctions[e]*ndim;
  int Nsize = nColSize*ndim;
  int NiSize = ndim*ndim;
  int NindexStart;
  double Ni;
  // Create N Matrix of size 3x24 in column major format.
  // N is computed for each Gauss Quadrature point
  double *N = (double*)calloc(Nsize, sizeof(double));
  for (int k = 0; k < GaussPoints[e]; k++) {
    // Populate N for each Gauss Point
    for (int n = 0; n < nShapeFunctions[e]; ++n) {
      NindexStart = n*NiSize;
      Ni = shp[gptr[e]+k*GaussPoints[e]+n];

      N[NindexStart] = Ni;
      N[NindexStart+4] = Ni;
      N[NindexStart+8] = Ni;
    }
    // Compute N^T N
    dgemm_(chy, chn, &nColSize, &nColSize, &ndim, &one, N, &ndim, \
        N, &ndim, &zero, MeGQ, &nColSize);
    // TODO(Anil) Gauss weights and det J product can be combined in shape
    // function
    int wIndex = gpPtr[e]+k;
    const double preFactor = gaussWeights[wIndex]*detJacobian[wIndex];
    // Me = \Sum_j w_j (N^T N Det(J))_j
    for (int n = 0; n < mLocalSize; ++n) {
      Me[n] += MeGQ[n]*preFactor;
    }
  }
  for (int n = 0; n < mLocalSize; ++n) {
    Me[n] *= rho;
  }
  // print Me Matrix
  if (debug && 1==0) {
    printf("DEBUG : Printing Me (Mass Matrix) for Element %d\n", e);
    for (int j = 0; j < bColSize; ++j) {
      for (int k = 0; k < bColSize; ++k) {
        printf("%.4f\t", Me[j+k*bColSize]);
      }
      printf("\n");
    }
  }
  free(MeGQ);
  free(N);
  return;
}
void LumpMassMatrix(void) {
  const int massSize = nnodes*ndim;
  // Lump mass matrix by summing up along the row
  for(int i = 1; i < massSize; ++i) {
    for(int j = 0; j < massSize; ++j) {
      mass[j] += mass[j+i*massSize];
    }
  }
	if(debug && 1==0){
	  printf("Lumped Mass\n");
	  for(int j = 0; j < massSize; ++j) {
	    printf("%.6f\n", mass[j]);
	  }
	}
}
