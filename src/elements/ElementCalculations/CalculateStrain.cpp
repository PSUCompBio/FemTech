#include "FemTech.h"
#include "blas.h"

#include <math.h>

const double kPi = 4.0*atan(1.0);

void CalculateMaximumPrincipalStrain(int elm, double* currentStrainMax, \
    double *currentStrainMin, double *currentShearMax) {
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

  double p1 = b*b + c*c + e*e;
  double eps1, eps2, eps3;
  if (p1 == 0) {
   // A is diagonal.
   eps1 = a;
   eps2 = d;
   eps3 = f;
  } else {
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

    eps1 = sqrtQ*cos(theta/3.0)+I1;
    eps2 = sqrtQ*cos((theta+2.0*kPi)/3.0)+I1;
    eps3 = sqrtQ*cos((theta+4.0*kPi)/3.0)+I1;
  }
  
  double min, max;
  max = fmax(eps3, fmax(eps2, eps1));
  min = fmin(eps3, fmin(eps2, eps1));
  *currentShearMax = 0.5*(max - min);
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
  double *E = &Eavg[elm*matSize];
    if (strcmp(ElementType[elm], "T3D2") == 0){
	double newlengthx = (coordinates[ndim*connectivity[eptr[elm]]+0]+displacements[ndim*connectivity[eptr[elm]]+0]-coordinates[ndim*connectivity[eptr[elm]+1]+0]-displacements[ndim*connectivity[eptr[elm]+1]+0]);
	double newlengthy = (coordinates[ndim*connectivity[eptr[elm]]+1]+displacements[ndim*connectivity[eptr[elm]]+1]-coordinates[ndim*connectivity[eptr[elm]+1]+1]-displacements[ndim*connectivity[eptr[elm]+1]+1]);
	double newlengthz = (coordinates[ndim*connectivity[eptr[elm]]+2]+displacements[ndim*connectivity[eptr[elm]]+2]-coordinates[ndim*connectivity[eptr[elm]+1]+2]-displacements[ndim*connectivity[eptr[elm]+1]+2]);
	double oldlengthx = (coordinates[ndim*connectivity[eptr[elm]]+0]-coordinates[ndim*connectivity[eptr[elm]+1]+0]);
	double oldlengthy = (coordinates[ndim*connectivity[eptr[elm]]+1]-coordinates[ndim*connectivity[eptr[elm]+1]+1]);
	double oldlengthz = (coordinates[ndim*connectivity[eptr[elm]]+2]-coordinates[ndim*connectivity[eptr[elm]+1]+2]);
	double newlength = newlengthx*newlengthx + newlengthy*newlengthy + newlengthz*newlengthz;
	double oldlength = oldlengthx*oldlengthx + oldlengthy*oldlengthy + oldlengthz*oldlengthz;
	E[0] = (newlength-oldlength)/2*oldlength;
	}
    else{  
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
}
