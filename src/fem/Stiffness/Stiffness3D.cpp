#include "FemTech.h"
#include "blas.h"

double *stiffness;

void AssembleStiffnessMatrix() {
	// set the debug flag for this file
	int debug = 1;
  // Create global stiffness matrix
  const int kSize = nnodes*ndim;
  stiffness = (double*)calloc(kSize*kSize, sizeof(double));

	int Csize = 6;

	// for debugging
	if (debug) {
		for (int i = 0; i < nelements; i++) {
			printf("shp array e.%d with %d shp functions\n", i, GaussPoints[i]);
			for (int k = 0; k < nShapeFunctions[i]; k++) {
				printf("int %d | dshp: %8.5f %8.5f %8.5f shp: %8.5f\n", k, dshp[dsptr[i] + k * ndim + 0],
														   dshp[dsptr[i] + k * ndim + 1],
														   dshp[dsptr[i] + k * ndim + 2],
														   shp[gptr[i] + k]);

			}
			printf("\n");
		}
	}

  for (int e = 0; e < nelements; ++e) {
    // number of shape functions * ndim
    int bColSize = nShapeFunctions[e]*ndim;
    int Bsize = bColSize*Csize;
    int keSize = bColSize*bColSize;
    // Create B Matrix of size 24x6 in column major format.
    // B is computed for each Gauss Quadrature point
    double *B = (double*)calloc(Bsize, sizeof(double));
    // Create local element stiffness matrix Ke of size 24x24 in column major
    // format
    double *Ke = (double*)calloc(keSize, sizeof(double));
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
    if (debug) {
      printf("RK:DEBUG : Printing Ke (Elemental Stiffness Matrix) for Element %d\n", e);
      for (int j = 0; j < bColSize; ++j) {
        for (int k = 0; k < bColSize; ++k) {
          printf("%.4f\t", Ke[j+k*bColSize]);
        }
        printf("\n");
      }
    }
    // Move local matrix to global matrix
    // TODO(Anil) remove the need for a local stiffness matrix by directly assigning
    // to global matrix
    for (int n = 0; n < nShapeFunctions[e]; ++n) {
      int g1Index = connectivity[eptr[e]+n];
      for (int m = 0; m < ndim; ++m) {
        for (int l = 0; l < nShapeFunctions[e]; ++l) {
          int g2Index = connectivity[eptr[e]+l];
          for (int k = 0; k < ndim; ++k) {
            stiffness[(g1Index*ndim+m)*kSize+g2Index*ndim+k] += Ke[(n*ndim+m)*bColSize+l*ndim+k];
          }
        }
      }
    }

    free(B);
    free(Ke);
    free(BtC);
    free(KeGQ);
  }
  if (debug) {
    printf("DEBUG : Printing Full Stiffness Matrix\n");
    for (int j = 0; j < kSize; ++j) {
      for (int k = 0; k < kSize; ++k) {
        printf("%.4f\t", stiffness[j+k*kSize]);
      }
      printf("\n");
    }
  }

	return;
}
