#include "FemTech.h"

/** For all elements -- this function calculates the minimum critical timestep */
double StableTimeStep() {

	double dtElem;
	double dtMin = huge;
  bool isNotRigid;
  // int minElementID;
  // int elemSkipped = 0;

	for(int i=0;i<nelements;i++){
    isNotRigid = false;
    for (int j = eptr[i]; j < eptr[i+1]; ++j) {
      int index = connectivity[j]*ndim;
      if (!(boundary[index]&boundary[index+1]&boundary[index+2])) {
        isNotRigid = true;
        break;
      }
    }
    if (isNotRigid) {
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

  if(dtMin < FailureTimeStep){
		FILE_LOG_MASTER(ERROR, "Timestep too small, dt = %15.9e", dtMin);
    TerminateFemTech(19);
	}
  return dtMin;
}
