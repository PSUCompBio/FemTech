#include "FemTech.h"

double *rhs;

void ApplySteadyBoundaryConditions(void) {
  //TODO(Anil) : Is elimination better ?
  // Currently the know prescribed displacement row and col are not eliminated
  // Prescribed displacement BC is incorporated into the equation
  const int matSize = nnodes*ndim;
  rhs = (double*)calloc(matSize, sizeof(double));
  for (int i = 0; i < nSpecifiedDispBC; ++i) {
    const int dof = nodeDOFDispBC[i];
    for (int j = 0; j < matSize; ++j) {
      stiffness[j*matSize+dof] = 0.0;
    }
    stiffness[matSize*dof+dof] = 1.0;
    rhs[dof] = nodeValueDispBC[i];
  }
  if (debug) {
    printf("DEBUG : Printing Full Stiffness Matrix\n");
    for (int j = 0; j < matSize; ++j) {
      for (int k = 0; k < matSize; ++k) {
        printf("%.4f\t", stiffness[j+k*matSize]);
      }
      printf("\n");
    }
    printf("DEBUG : Printing RHS\n");
    for (int j = 0; j < matSize; ++j) {
      printf("%.4f\n", rhs[j]);
    }
  }
}
