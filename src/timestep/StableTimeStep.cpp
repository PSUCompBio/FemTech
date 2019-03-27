#include "FemTech.h"
#include "blas.h"

/** For all elements -- this function calculates the minimum critical timestep */
double StableTimeStep() {

	double dt;
	double dtMin = HUGE;

	for(int i=0;i<nelements;i++){
 		dt = CalculateTimeStep(i);
		if (dt < dtMin){
			dtMin = dt;
		}
	}

  if(dtMin < FailureTimeStep){
		printf("Simulation Failed - Timestep too small.\n");
		printf("Timestep is: %3.3e\n", dt);
	}

  return dtMin;
}
