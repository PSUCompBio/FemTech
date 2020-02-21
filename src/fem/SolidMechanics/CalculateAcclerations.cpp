#include "FemTech.h"
#include "blas.h"

void CalculateAccelerations(){
  FILE_LOGArray(DEBUGLOGIGNORE, f_net, nDOF, "RHS of Solution");

  for (int i = 0; i < nDOF; ++i) {
    if (!boundary[i]) {
      accelerations[i] = f_net[i]/mass[i];
    }
  }
	return;
}
