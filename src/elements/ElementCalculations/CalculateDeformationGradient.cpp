#include "FemTech.h"
#include "blas.h"

void CalculateDeformationGradient(int e, int gp) {
  // following Bonet and Wood; F = xlocal*DN_X in flagshyp
  if ((eptr[e + 1] - eptr[e]) != nShapeFunctions[e]) {
    FILE_LOG_SINGLE(WARNING, "CalculateDeformationGradient.cpp: (eptr[e+1]-eptr[e]) "
           "!= nShapeFunctions[e])\nCheck it out, probally a bug in allocation");
  }
  double theSum = 0.0;
  for (int i = 0; i < ndim; i++) {
    for (int j = 0; j < ndim; j++) {
      theSum = 0.0;
      int index = fptr[e] + ndim * ndim * gp + ndim * j + i;
      for (int k = 0; k < nShapeFunctions[e]; k++) {
        int a = eptr[e];
        int node_a = connectivity[a + k];
        int index2 = dsptr[e] + gp * GaussPoints[e] * ndim + k * ndim + j;
        theSum = theSum + (coordinates[ndim * node_a + i] +
                           displacements[ndim * node_a + i]) *
                              dshp[index2];
      } // loop on k
      F[index] = theSum;
    }
  }

  FILE_LOGMatrix_SINGLE(DEBUGLOGIGNORE, &F[fptr[e] + ndim * ndim * gp], ndim, \
      ndim, "--------F for gauss point %d --------", gp);
  return;
}
