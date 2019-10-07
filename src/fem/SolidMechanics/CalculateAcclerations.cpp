#include "FemTech.h"
#include "blas.h"

void CalculateAccelerations(){
  FILE_LOGArray(DEBUGLOGIGNORE, f_net, ndim*nnodes, "RHS of Solution");

  for (int i = 0; i < ndim*nnodes; ++i) {
    if (!boundary[i]) {
      accelerations[i] = f_net[i]/mass[i];
    }
  }
	return;
}
