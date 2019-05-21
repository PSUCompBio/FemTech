#include "FemTech.h"
#include "blas.h"

void CalculateDeformationGradient(int e, int gp) {

  // good example of how to reference deformation gradient, F
  if (debug && 1 == 1) {
    printf("shp array e.%d with %d Gauss points, each with %d shp functions \n",
           e, GaussPoints[e], nShapeFunctions[e]);
    for (int k = 0; k < nShapeFunctions[e]; k++) {
      printf(" shp: %4.4f dshp: %8.4f %8.4f %8.4f\n",
             shp[gptr[e] + gp * GaussPoints[e] + k],
             dshp[dsptr[e] + gp * GaussPoints[e] * ndim + k * ndim + 0],
             dshp[dsptr[e] + gp * GaussPoints[e] * ndim + k * ndim + 1],
             dshp[dsptr[e] + gp * GaussPoints[e] * ndim + k * ndim + 2]);
    }
  }

  // following Bonet and Wood; F = xlocal*DN_X in flagshyp
  if ((eptr[e + 1] - eptr[e]) != nShapeFunctions[e]) {
    printf("Warning from CalculateDeformationGradient.cpp: (eptr[e+1]-eptr[e]) "
           "!= nShapeFunctions[e])");
    printf("Check it out, probally a bug in allocation");
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
        if (i==0 && j ==0) {
          printf("F11 : coord : %12.8f, displ : %12.8f, total : %12.8f\n", coordinates[ndim * node_a + i], displacements[ndim * node_a + i], (coordinates[ndim * node_a + i] + displacements[ndim * node_a + i]));
        }
        if (i==2 && j ==2) {
          printf("F33 : coord : %12.8f, displ : %12.8f, total : %12.8f\n", coordinates[ndim * node_a + i], displacements[ndim * node_a + i], (coordinates[ndim * node_a + i] + displacements[ndim * node_a + i]));
        }
      } // loop on k
      F[index] = theSum;
    }
  }

  if (debug && 1 == 1) {
    printf("--------F for gauss point %d --------\n", gp);
    for (int i = 0; i < ndim; i++) {
      for (int j = 0; j < ndim; j++) {
        int index = fptr[e] + ndim * ndim * gp + ndim * j + i;
        printf("%3.3e   ", F[index]);
      } // loop on j
      printf("\n");
    } // loop on i
  }   // if debug

  return;
}
