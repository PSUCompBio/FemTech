#include "FemTech.h"
#include "blas.h"

void CalculateAccelerations(){
  FILE_LOGArray(DEBUGLOGIGNORE, f_net, nDOF, "RHS of Solution");
  for (int i = 0; i < nDOF; ++i) {
    accelerations[i] = (f_net[i]-dynamicDamping*velocities_half[i])/mass[i];
  }
	return;
}


