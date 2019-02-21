#include "FemTech.h"

#include "blas.h"

double* displacement;

void SolveSteadyImplicit(void) {
  int matSize = nnodes*ndim-nSpecifiedDispBC;
  displacement = (double*)malloc(matSize*sizeof(double));
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
  memcpy(displacement, rhs, matSize*sizeof(double));
  if (debug) {
    printf("DEBUG : Printing Displacement Solution\n");
    for (int j = 0; j < matSize; ++j) {
      printf("%12.4f\n", displacement[j]);
    }
  }
  return;
}
