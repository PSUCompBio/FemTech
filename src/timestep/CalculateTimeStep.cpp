#include "FemTech.h"
#include "blas.h"

/** For a single element - this function calculates the volume of the element
    and calculates the critical time step based on the wave speed.*/

double CalculateTimeStep(int e) {
	double dtElem;
 	// characteristic element length
  double le = CalculateCharacteristicLength(e);
	//wave speed of material
  int pide = pid[e];
	double ce = waveSpeed[pide];
	dtElem = le/ce;
  // FILE_LOG_MASTER(WARNING, "Ce : %3.3e, Le : %3.3e, dt = %3.3e", ce, le, dtElem);
	return dtElem;
}
