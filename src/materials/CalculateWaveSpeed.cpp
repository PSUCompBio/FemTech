#include "FemTech.h"

#include <math.h>

double CalculateWaveSpeed(const unsigned int partID) {
  double E = 0.0;
  double *partProperties = &(properties[partID*MAXMATPARAMS]);
  const double rho = partProperties[0];
  const unsigned int matID = materialID[partID];

  // Variables required inside switch
  double lambda, mu, K, muI, viscoEffect;
  unsigned int nTerm, nPronyL;
  switch (matID) {
    case 0 : // Rigid-body motion
             break;
    case 1 : //Compressible Neohookean
    case 2 : // St. Venant-Kirchhoff
    case 3 : // Linear Elastic
    case 4 : // HGO with isotropic fiber distribution
             lambda = partProperties[2];
             mu = partProperties[1];
             break;
    case 5 : // HGO with isotropic fiber distribution and viscoelasticity
             lambda = partProperties[2];
             muI = partProperties[1];
             nPronyL = static_cast<int>(partProperties[5]);
             viscoEffect = 1.0;
             for (int i = 0; i < nPronyL; ++i) {
               const double gi = partProperties[6+2*i];
               const double taui = partProperties[7+2*i];
               viscoEffect = viscoEffect + gi*exp(-1.0*(Time-tInitial)/taui);
             }
             mu = muI * viscoEffect;
             break;
    case 6 : // LS-Dyna Maxwell viscoelastic material equivalent
             lambda = partProperties[1];
             muI = partProperties[2];
             viscoEffect = 1.0 + partProperties[3]*exp(-1.0*(Time-tInitial)*properties[4]);
             mu = muI *viscoEffect;
             break;
    case 7 : // Ogden
             K = partProperties[1];
             nTerm = static_cast<int>(partProperties[2]);
             mu = 0.0;
             for (int i = 0; i < nTerm; ++i) {
               mu = mu + partProperties[3+2*i]*partProperties[4+2*i];
             }
             mu = 0.5*mu;
             lambda = K - 2.0*mu/3.0;
             break;
    case 8 : // Ogden Viscoelastic
             K = partProperties[1];
             // Compute mu_\infty = 0.5*\Sigma alpha_i mu_i
             nTerm = static_cast<int>(partProperties[2]);
             muI = 0.0;
             for (int i = 0; i < nTerm; ++i) {
               muI = muI + partProperties[3+2*i]*partProperties[4+2*i];
             }
             muI = 0.5*muI;
             // Compute mu = mu_\infity*(1+\Sigma g_i*exp(-t/tau_i)
             nPronyL = static_cast<int>(partProperties[3+2*nTerm]);
             viscoEffect = 1.0;
             for (int i = 0; i < nPronyL; ++i) {
               const double gi = partProperties[4+2*(i+nTerm)];
               const double taui = partProperties[5+2*(i+nTerm)];
               viscoEffect = viscoEffect + gi*exp(-1.0*(Time-tInitial)/taui);
             }
             mu = muI * viscoEffect;
             lambda = K - 2.0*mu/3.0;
             break;
    default : FILE_LOG_SINGLE(ERROR, "Unknown material type in wave speed computation");
              TerminateFemTech(1);
  }
  // Compute dilatational wave speed
  const double ce = sqrt((lambda+2.0*mu)/rho);
  // FILE_LOG_SINGLE(WARNING, "Part : %d, Material type, %d with wave speed %15.9f", partID, matID, ce);
  return ce;
}
