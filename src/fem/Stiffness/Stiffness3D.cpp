#include "FemTech.h"
#include "blas.h"

void StiffnessElementMatrix(double* Ke, int e) {
	int Csize = 6;
  // number of shape functions * ndim
  int bColSize = nShapeFunctions[e]*ndim;
  int Bsize = bColSize*Csize;
  int keSize = bColSize*bColSize;
  // Create B Matrix of size 24x6 in column major format.
  // B is computed for each Gauss Quadrature point
  for (int k = 0; k < Bsize; ++k) {
    B[k] = 0.0;
  }
  // Create temperory variables to store intermediate results
  double *BtC = (double*)calloc(Bsize, sizeof(double));
  double *KeGQ = (double*)calloc(keSize, sizeof(double));

  // Computing Ke matrix
  int BiSize = ndim*Csize;
  double dNdx, dNdy, dNdz;
  int indexStart, BindexStart;
  for (int k = 0; k < GaussPoints[e]; k++) {
    // Populate B for each Gauss Point
    for (int n = 0; n < nShapeFunctions[e]; ++n) {
      indexStart = dsptr[e]+(k*GaussPoints[e]+n)*ndim;
      BindexStart = n*BiSize;
      dNdx = dshp[indexStart];
      dNdy = dshp[indexStart+1];
      dNdz = dshp[indexStart+2];

      B[BindexStart] = dNdx;
      B[BindexStart+7] = dNdy;
      B[BindexStart+14] = dNdz;

      B[BindexStart+9] = dNdz;
      B[BindexStart+15] = dNdy;
      B[BindexStart+4] = dNdz;
      B[BindexStart+16] = dNdx;
      B[BindexStart+5] = dNdy;
      B[BindexStart+11] = dNdx;
    }
    // Compute B^T C
    dgemm_(chy, chn, &bColSize, &Csize, &Csize, &one, B, &Csize, \
        C, &Csize, &zero, BtC, &bColSize);
    // Compute B^T C B
    dgemm_(chn, chn, &bColSize, &bColSize, &Csize, &one, BtC, &bColSize, \
        B, &Csize, &zero, KeGQ, &bColSize);
    int wIndex = gpPtr[e]+k;
    const double preFactor = gaussWeights[wIndex]*detJacobian[wIndex];
    // Ke = \Sum_j w_j (B^T C B Det(J))_j
    for (int n = 0; n < keSize; ++n) {
      Ke[n] += KeGQ[n]*preFactor;
    }
  }

  // print Ke Matrix
#ifdef DEBUG
  if (debug && 1==0) {
    printf("RK:DEBUG : Printing Ke (Elemental Stiffness Matrix) for Element %d\n", e);
    for (int j = 0; j < bColSize; ++j) {
      for (int k = 0; k < bColSize; ++k) {
        printf("%12.4f", Ke[j+k*bColSize]);
      }
      printf("\n");
    }
  }
#endif //DEBUG
  free(B);
  free(BtC);
  free(KeGQ);
}
