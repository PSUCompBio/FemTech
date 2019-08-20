#include "FemTech.h"
#include "blas.h"

/** For a single element - this function calculates the volume of the element
    and calculates the critical time step based on the wave speed.*/

double CalculateTimeStep(int e) {
	double dt;
 	// characteristic element length
  double le = CalculateCharacteristicLength(e);
	//wave speed of material
  int pide = pid[e];
  double mu = properties[MAXMATPARAMS * pide + 1];
  double lambda = properties[MAXMATPARAMS * pide + 2];
  double rho = properties[MAXMATPARAMS * pide + 0];
  double nu = 0.5*lambda/(lambda+mu);
  // ce calculated using dilatatoin wave speed
	double ce = sqrt(lambda*(1.0/nu-1.0)/rho);
	dt = le/ce;
	return dt;
}
