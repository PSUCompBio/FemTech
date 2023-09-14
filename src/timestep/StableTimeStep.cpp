#include "FemTech.h"

/** For all elements -- this function calculates the minimum critical timestep */
double StableTimeStep() {
	double dtElem;
	double dtMin = huge;
  bool isNotRigid;
  // int minElementID;
  // int elemSkipped = 0;
  
  // Update wave speed of viscoelastic parts
  for (int i = 0; i < nPIDglobal; ++i) {
    if(viscoElasticPart[i]) {
      waveSpeed[i] = CalculateWaveSpeed(i);
    }
  }

	for (int i = 0; i < nelements; ++i) {
    if (materialID[pid[i]] != 0)
	if(pid[i] != 10) {
      dtElem = CalculateTimeStep(i);
      if (dtElem < dtMin){
        dtMin = dtElem;
        // minElementID = i;
      }
    } 
    // else {
    //   elemSkipped = elemSkipped + 1;
    // }
	}
  // FILE_LOG(DEBUGLOG, "Minimum dt from element %d", minElementID);
  // FILE_LOG(DEBUGLOG, "Elements skipped %d", elemSkipped);
  MPI_Allreduce(MPI_IN_PLACE, &dtMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  if(dtMin < FailureTimeStep) {
		FILE_LOG_MASTER(ERROR, "Timestep too small, dt = %15.9e", dtMin);
    // Find the element causing low time step and report the element
    // on log before termination
    for (int i = 0; i < nelements; ++i) {
      // materialID = 0 => Rigid element
      if (materialID[pid[i]] != 0) 
	if(pid[i] != 10) {
        dtElem = CalculateTimeStep(i);
        if (dtElem < FailureTimeStep) {
		      FILE_LOG_SINGLE(ERROR, "Small timestep detected in element number : %8d, dt = %15.9e, wave speed = %15.9e", global_eid[i], dtElem, waveSpeed[i]);
        }
      } 
    }
	}
  if(dtMin > MaxTimeStep){
    dtMin = MaxTimeStep;
	}
  return dtMin;
}
