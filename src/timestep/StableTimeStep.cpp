#include "FemTech.h"
#include "blas.h"

/** For all elements -- this function calculates the minimum critical timestep */
double StableTimeStep() {

	double dt;
	double dtMin = huge;
  bool isNotRigid;
  // int minElementID;
  // int elemSkipped = 0;

	for(int i=0;i<nelements;i++){
    isNotRigid = false;
    for (int j = eptr[i]; j < eptr[i+1]; ++j) {
      if (!boundary[connectivity[j]]) {
        isNotRigid = true;
        break;
      }
    }
    if (isNotRigid) {
      dt = CalculateTimeStep(i);
      if (dt < dtMin){
        dtMin = dt;
        // minElementID = i;
      }
    } 
    // else {
    //   elemSkipped = elemSkipped + 1;
    // }
	}
  // printf("INFO(%d): Minimum dt from element %d\n", world_rank, minElementID);
  // printf("INFO(%d): Elements skipped %d\n", world_rank, elemSkipped);
  MPI_Allreduce(MPI_IN_PLACE, &dtMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  if(dtMin < FailureTimeStep){
		printf("Simulation Failed - Timestep too small.\n");
		printf("Timestep is: %3.3e\n", dtMin);
	}

  return dtMin;
}
