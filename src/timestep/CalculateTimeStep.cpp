#include "FemTech.h"
#include "blas.h"

/** For a single element - this function calculates the volume of the element
    and calculates the critical time step based on the wave speed.*/

double CalculateTimeStep(int e) {
	double dt;
  const double dtMax = 0.1;
 	// characteristic element length
  double le = 1.0;
	//wave speed of material
  double mu = properties[MAXMATPARAMS * e + 1];
  double lambda = properties[MAXMATPARAMS * e + 2];
  double rho = properties[MAXMATPARAMS * e + 0];
  double nu = 0.5*lambda/(lambda+mu);
  // ce calculated using dilatatoin wave speed
	double ce = sqrt(lambda*(1.0/nu-1.0)/rho);
	dt = le/ce;
  if (dt > dtMax) {
    dt = dtMax;
  }
  printf("Time Step : %.24f\n", dt);
	return dt;
}
