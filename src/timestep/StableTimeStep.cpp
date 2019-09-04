#include "FemTech.h"
#include "blas.h"

/** For all elements -- this function calculates the minimum critical timestep */
double StableTimeStep() {

	double dt;
	double dtMin = huge;

	for(int i=0;i<nelements;i++){
 		dt = CalculateTimeStep(i);
		if (dt < dtMin){
			dtMin = dt;
		}
	}
  MPI_Allreduce(MPI_IN_PLACE, &dtMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  if(dtMin < FailureTimeStep){
		printf("Simulation Failed - Timestep too small.\n");
		printf("Timestep is: %3.3e\n", dtMin);
	}

  return dtMin;
}
