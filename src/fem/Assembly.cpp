#include "FemTech.h"

double *stiffness;
double *mass;

void Assembly(char *operation) {
	if (strcmp("mass", operation) == 0) {
    // Create global mass matrix
    const int massSize = nnodes*ndim;
    mass = (double*)calloc(massSize*massSize, sizeof(double));
    for (int e = 0; e < nelements; ++e) {
      int bColSize = nShapeFunctions[e]*ndim;
      int mLocalSize = bColSize*bColSize;
      double *Me = (double*)calloc(mLocalSize, sizeof(double));
      MassElementMatrix(Me, e);
      // Move local matrix to global matrix
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
    }
#ifdef DEBUG
    if (debug && 1==0) {
      printf("DEBUG : Printing Full Mass Matrix\n");
      for (int j = 0; j < massSize; ++j) {
        for (int k = 0; k < massSize; ++k) {
          printf("%.4f\t", mass[j+k*massSize]);
        }
        printf("\n");
      }
    }
#endif //DEBUG
	} else {
	  if (strcmp("stiffness", operation) == 0) {
      // Create global stiffness matrix
      const int kSize = nnodes*ndim;
      stiffness = (double*)calloc(kSize*kSize, sizeof(double));
      for (int e = 0; e < nelements; ++e) {
        // number of shape functions * ndim
        int bColSize = nShapeFunctions[e]*ndim;
        int keSize = bColSize*bColSize;
        // Create local element stiffness matrix Ke of size 24x24 in column major
        // format
        double *Ke = (double*)calloc(keSize, sizeof(double));
        StiffnessElementMatrix(Ke, e);
        // Move local matrix to global matrix
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
        free(Ke);
      }
#ifdef DEBUG
      if (debug && 1==0) {
        printf("DEBUG : Printing Full Stiffness Matrix\n");
        for (int j = 0; j < kSize; ++j) {
          for (int k = 0; k < kSize; ++k) {
            printf("%12.4f", stiffness[j+k*kSize]);
          }
          printf("\n");
        }
      }
#endif //DEBUG
    }
  }
	return;
}
