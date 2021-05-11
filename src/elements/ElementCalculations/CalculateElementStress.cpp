#include "FemTech.h"

void CalculateElementStress(unsigned int e, double* stress) {
  unsigned int limit = ndim*(ndim+1)/2;
  for (int i = 0; i < limit; ++i) {
    stress[i] = 0.0;
  }
  const int countGP = GaussPoints[e];
  double gwSum = 0.0;
  for(int gp = 0; gp < countGP; ++gp) {
    double * sigmaLocal = &(sigma_n[sigmaptr[e]+6*gp]);
    const unsigned int wIndex = gpPtr[e]+gp;
    const double gW = gaussWeights[wIndex];
    gwSum = gwSum + gW;
    for (int i = 0; i < limit; ++i) {
      stress[i] = stress[i] + gW*sigmaLocal[i];
    }
  }
  for (int i = 0; i < limit; ++i) {
    stress[i] = stress[i]/gwSum;
  }
}
