#include "FemTech.h"

#include "blas.h"

double *mass;

void AssembleMassMatrix() {
  // Create global mass matrix
  const int massSize = nnodes*ndim;
  mass = (double*)calloc(massSize*massSize, sizeof(double));

	int Csize = 6;

  for (int e = 0; e < nelements; ++e) {
    // number of shape functions * ndim
    int bColSize = nShapeFunctions[e]*ndim;
    int Bsize = bColSize*Csize;
    int mLocalSize = bColSize*bColSize;

    //TODO(Anil) remove allocations inside the loop
	  double *Me = (double*)calloc(mLocalSize, sizeof(double));
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
    if (debug) {
      printf("DEBUG : Printing Me (Mass Matrix) for Element %d\n", e);
      for (int j = 0; j < bColSize; ++j) {
        for (int k = 0; k < bColSize; ++k) {
          printf("%.4f\t", Me[j+k*bColSize]);
        }
        printf("\n");
      }
    }
    // Move local matrix to global matrix
    // TODO(Anil) remove the need for a local mass matrix by directly assigning
    // to global matrix
    for (int n = 0; n < nShapeFunctions[e]; ++n) {
      int g1Index = connectivity[eptr[e]+n];
      for (int m = 0; m < ndim; ++m) {
        for (int l = 0; l < nShapeFunctions[e]; ++l) {
          int g2Index = connectivity[eptr[e]+l];
          for (int k = 0; k < ndim; ++k) {
            mass[(g1Index*ndim+m)*massSize+g2Index*ndim+k] += Me[(n*ndim+m)*bColSize+l*ndim+k];
          }
        }
      }
    }
    free(Me);
    free(MeGQ);
    free(N);
  }
  printf("Mass size : %d\n", massSize);
  if (debug) {
    printf("DEBUG : Printing Full Mass Matrix\n");
    for (int j = 0; j < massSize; ++j) {
      for (int k = 0; k < massSize; ++k) {
        printf("%.4f\t", mass[j+k*massSize]);
      }
      printf("\n");
    }
  }

	return;
}
