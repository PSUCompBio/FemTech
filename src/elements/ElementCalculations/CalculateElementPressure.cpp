#include "FemTech.h"

double CalculateElementPressure(unsigned int e) {
  double pressure = 0.0;
  const int countGP = GaussPoints[e];
  double gwSum = 0.0;
  for(int gp = 0; gp < countGP; ++gp) {
    double * sigmaLocal = &(sigma_n[sigmaptr[e]+6*gp]);
    const unsigned int wIndex = gpPtr[e]+gp;
    const double gW = gaussWeights[wIndex];
    gwSum = gwSum + gW;
    for (int i = 0; i < ndim; ++i) {
      pressure = pressure + gW*sigmaLocal[i];
    }
  }
  pressure = -pressure/3.0/gwSum;
  return pressure;
}
