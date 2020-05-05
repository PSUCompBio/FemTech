#include "FemTech.h"
#include "blas.h"

#include <math.h>

const double kPi = 4.0*atan(1.0);

void CalculateMaximumPrincipalStrain(int elm, double* currentStrainMax, \
    double *currentStrainMin) {
  const int matSize = ndim*ndim;

  double *E = &Eavg[elm*matSize];
  for (int i = 0; i < matSize; ++i) {
    E[i] = 0.0;
  }
  const int countGP = GaussPoints[elm];
  double preFactor = 0.5/((double)countGP);
  for(int gp = 0; gp < countGP; ++gp) {
    int index = fptr[elm] + matSize * gp;
    double *F_element_gp = &(F[index]);
    // Compute Green-Lagrange Tensor: E= (1/2)*(F^T*F - I)
    dgemm_(chy, chn, &ndim, &ndim, &ndim, &preFactor, F_element_gp, &ndim,
           F_element_gp, &ndim, &one, E, &ndim);
  }
  E[0] -= 0.5;
  E[4] -= 0.5;
  E[8] -= 0.5;

  // Calculate the eigen values
  // E is symm
  double a, b, c, d, e, f;

	a = E[0]; b = E[1]; c = E[2];
	d = E[4]; e = E[5]; f = E[8];

	// charcteristic eq. : x^3 - I1*x^2 + I2*x -I3 = 0
	double I1 = a + d + f;
	double I2 = a*(d+f)+d*f-b*b-c*c-e*e;
	double I3 = a*d*f+2.0*b*c*e-b*b*f-c*c*d-e*e*a;

  // https://www.continuummechanics.org/principalstress.html
  double Q = (3.0*I2-I1*I1)/9.0;
  double R = (2.0*I1*I1*I1-9.0*I1*I2+27.0*I3)/54.0;
  double theta = acos(R/sqrt(-Q*Q*Q));
  double sqrtQ = 2.0*sqrt(-Q);
  I1 = I1/3.0;

  double eps1, eps2, eps3;
  eps1 = sqrtQ*cos(theta/3.0)+I1;
  eps2 = sqrtQ*cos((theta+2.0*kPi)/3.0)+I1;
  eps3 = sqrtQ*cos((theta+4.0*kPi)/3.0)+I1;

  double min, max;
  max = fmax(eps3, fmax(eps2, eps1));
  min = fmin(eps3, fmin(eps2, eps1));
  if (max > 0.0) {
    *currentStrainMax = max;
  } else {
    *currentStrainMax = 0.0;
  }
  if (min < 0.0) {
    *currentStrainMin = min;
  } else {
    *currentStrainMin = 0.0;
  }
}

void CalculateStrain() {
  const int matSize = ndim*ndim;
  for (int i = 0; i < matSize*nelements; ++i) {
    Eavg[i] = 0.0;
  }
  for (int elm = 0; elm < nelements; ++elm) {
    /*The next line says local E is part of Eavg and starts from elm*matSize.
    So when we assign to E, we are in effect assigning to Eavg.*/
    double *E = &Eavg[elm*matSize];
    const int countGP = GaussPoints[elm];
    double preFactor = 0.5/((double)countGP);
    for(int gp = 0; gp < countGP; ++gp) {
      int index = fptr[elm] + matSize * gp;
      double *F_element_gp = &(F[index]);
      // Compute Green-Lagrange Tensor: E= (1/2)*(F^T*F - I)
      dgemm_(chy, chn, &ndim, &ndim, &ndim, &preFactor, F_element_gp, &ndim,
            F_element_gp, &ndim, &one, E, &ndim);
    }
    E[0] -= 0.5;
    E[4] -= 0.5;
    E[8] -= 0.5;
  }
}
