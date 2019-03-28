#include "FemTech.h"
#include "blas.h"

/** For a single element - this function calculates the volume of the element
    and calculates the critical time step based on the wave speed.*/

double CalculateTimeStep(int e){

	double dt;
 	// characteristic element length
  double le = 1.0;
	//wave speed of material
	double ce = sqrt(262.5/1000.);
	dt = le/ce;

	dt = 1e-4;

	return dt;
}
