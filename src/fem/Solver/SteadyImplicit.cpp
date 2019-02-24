#include "FemTech.h"

#include "blas.h"

void SolveSteadyImplicit(void) {
  int matSize = nnodes*ndim-nSpecifiedDispBC;
  // To store pivot order in LU
  int* pivot = (int*)malloc(matSize*sizeof(int));
  int info;
  int oneI = 1;
  // Solve K u = f directly by using LU Decomposition
  dgesv_(&matSize, &oneI, stiffness, &matSize, pivot, rhs, &matSize, &info);
  if (info) {
    printf("LU Decomposition and solution failed with info code %d\n", info);
  }
  // Copy solution to displacement vector
  int j = 0;
  for (int i = 0; i < nnodes*ndim; ++i) {
    if(boundary[i] == 0) {
      displacements[i] = rhs[j];
      j = j + 1;
    }
  }

  if (debug) {
    printf("DEBUG : Printing Displacement Solution\n");
    for (int i = 0; i < nnodes; ++i) {
      for (int j = 0; j < ndim; ++j) {
        printf("%12.4f", displacements[i*ndim+j]);
      }
      printf("\n");
    }
  }
  return;
}
