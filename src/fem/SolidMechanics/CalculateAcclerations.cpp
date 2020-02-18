#include "FemTech.h"
#include "blas.h"

void CalculateAccelerations(){
#ifdef DEBUG
  if (debug && 1 == 0) {
    printf("RHS of solution\n");
    for (int i = 0; i < nDOF; ++i) {
	    printf("%d  %12.6f\n", i, f_net[i]);
	  }
  }
#endif //DEBUG
  for (int i = 0; i < nDOF; ++i) {
    if (!boundary[i]) {
      accelerations[i] = f_net[i]/mass[i];
    }
  }
	return;
}
