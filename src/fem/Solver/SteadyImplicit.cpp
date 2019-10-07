#include "FemTech.h"
#include "blas.h"

void SolveSteadyImplicit(void) {
  int matSize = nnodes*ndim-nSpecifiedDispBC;
  // To store pivot order in LU
  int* pivot = (int*)malloc(matSize*sizeof(int));
  int info;
  // Solve K u = f directly by using LU Decomposition
  dgesv_(&matSize, &oneI, stiffness, &matSize, pivot, rhs, &matSize, &info);
  if (info) {
    FILE_LOG_SINGLE(ERROR, "LU Decomposition and solution failed with info code %d", info);
  }
  // Copy solution to displacement vector
  int j = 0;
  for (int i = 0; i < nnodes*ndim; ++i) {
    if(boundary[i] == 0) {
      displacements[i] = rhs[j];
      j = j + 1;
    }
  }
  FILE_LOGMatrixRM(DEBUGLOGIGNORE, displacements, nnodes, ndim, \
      "Printing Displacement Solution");
  free(pivot);
  return;
}
