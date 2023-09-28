#include "FemTech.h"
#include "blas.h"

/** For a single element - this function calculates the volume of the element
    and calculates the critical time step based on the wave speed.*/

double CalculateTimeStep(int e) {
	double dtElem;
	double le;
  int pide = pid[e];
	double ce = waveSpeed[pide]; //wavespeed
 	// characteristic element length
	if(embed){
	  if (strcmp(ElementType[e], "T3D2") != 0){
	    le = CalculateCharacteristicLength(e);
	    dtElem = le/ce;
	  }
	  else dtElem = 10; //to ensure that the timestep is not affected by fibers
	} else{
  	  le = CalculateCharacteristicLength(e);
	  dtElem = le/ce;
	}
  // FILE_LOG_MASTER(WARNING, "Ce : %3.3e, Le : %3.3e, dt = %3.3e", ce, le, dtElem);
	return dtElem;
}
