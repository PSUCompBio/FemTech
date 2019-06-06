#include "FemTech.h"
#include "blas.h"

void CalculateAccelerations(){
	//printf("Calculate Accelerations.\n");
	//Accelerations = VectorXd::Zero(F_net.size());

  if (debug && 1 == 0) {
    printf("RHS of solution\n");
    for (int i = 0; i < ndim*nnodes; ++i) {
	    printf("%d  %12.6f\n", i, f_net[i]);
	  }
  }
  for (int i = 0; i < ndim*nnodes; ++i) {
    accelerations[i] = f_net[i]/mass[i];
    // printf("Acceleration %d : %.5f\n", i, accelerations[i]);
  }
  for (int i = 0; i < nnodes; ++i) {
    printf("Z Acceleration %d : %12.6f, %12.6f, %12.6f\n", i, accelerations[3*i+2], f_net[3*i+2], mass[3*i+2]);
  }
  for (int i = 0; i < nnodes; ++i) {
    printf("X Acceleration %d : %12.6f, %12.6f, %12.6f\n", i, accelerations[3*i], f_net[3*i], mass[3*i]);
  }
	return;
}
