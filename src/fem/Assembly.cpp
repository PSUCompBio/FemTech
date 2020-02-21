#include "FemTech.h"

double *stiffness;
double *mass;

void Assembly(char *operation) {
	if (strcmp("mass", operation) == 0) {
    // Create global mass matrix
    mass = (double*)calloc(nDOF*nDOF, sizeof(double));
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
              mass[(g1Index*ndim+m)*nDOF+g2Index*ndim+k] += Me[(n*ndim+m)*bColSize+l*ndim+k];
            }
          }
        }
      }
      free(Me);
    }
    FILE_LOGMatrix_SINGLE(DEBUGLOG, mass, nDOF, nDOF, \
        "Printing Full Mass Matrix");
	} else {
	  if (strcmp("stiffness", operation) == 0) {
      // Create global stiffness matrix
      stiffness = (double*)calloc(nDOF*nDOF, sizeof(double));
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
                stiffness[(g1Index*ndim+m)*nDOF+g2Index*ndim+k] += Ke[(n*ndim+m)*bColSize+l*ndim+k];
              }
            }
          }
        }
        free(Ke);
      }
      FILE_LOGMatrix_SINGLE(DEBUGLOG, stiffness, nDOF, nDOF, \
          "Printing Full Stiffness Matrix");
    }
  }
	return;
}
