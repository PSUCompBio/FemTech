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
  double rho = properties[MAXMATPARAMS * pide + 0];
  double mu, lambda;
  unsigned int matID = materialID[pide];
  // TODO precompute material wave speed and store for all moterial types used
  // in the simulation.
  if (matID == 6) { // Viscoelastic
    mu = properties[MAXMATPARAMS * pide + 2];
    double K = properties[MAXMATPARAMS * pide + 1];
    lambda = K + 2.0*mu/3.0;
  } else {
    if (matID == 7 || matID == 8) { // Ogden model
      // TODO : mu_1 used, not ideal. 
      mu = properties[MAXMATPARAMS * pide + 4];
      double K = properties[MAXMATPARAMS * pide + 1];
      lambda = mu + 2.0*K/3.0;
    } else {
      mu = properties[MAXMATPARAMS * pide + 1];
      lambda = properties[MAXMATPARAMS * pide + 2];
    }
  }
  FILE_LOG_SINGLE(DEBUGLOG, "Element %d, pid %d, Mu = %f, Lambda = %f", e, pide, mu, lambda);
  double nu = 0.5*lambda/(lambda+mu);
  // ce calculated using dilatatoin wave speed
	double ce = sqrt(lambda*(1.0/nu-1.0)/rho);
	dtElem = le/ce;
	return dtElem;
}
