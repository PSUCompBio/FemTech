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
  double mu = properties[MAXMATPARAMS * pide + 1];
  double lambda = properties[MAXMATPARAMS * pide + 2];
  double rho = properties[MAXMATPARAMS * pide + 0];
  double nu = 0.5*lambda/(lambda+mu);
  // ce calculated using dilatatoin wave speed
  // TODO : Store Element Ce during pre-processing
	double ce = sqrt(lambda*(1.0/nu-1.0)/rho);

	dtElem = le/ce;
  // FILE_LOG_MASTER(WARNING, "Ce : %3.3e, Le : %3.3e, dt = %3.3e", ce, le, dtElem);
	return dtElem;
}
