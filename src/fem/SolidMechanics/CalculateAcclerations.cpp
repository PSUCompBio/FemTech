#include "FemTech.h"
#include "blas.h"

void CalculateAccelerations(){
	//printf("Calculate Accelerations.\n");
	//Accelerations = VectorXd::Zero(F_net.size());

	//for (int i = 0; i < F_net.size(); i++) {
	//	Accelerations(i) = (F_net(i) / mm(i));
	//}
  for (int i = 0; i < ndim*nnodes; ++i) {
    accelerations[i] = f_net[i]/mass[i];
    printf("Acceleration %d : %.5f\n", i, accelerations[i]);
  }
	return;
}
