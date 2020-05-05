#include "FemTech.h"
#include "blas.h"

#include <math.h>

const double kPi = 4.0*atan(1.0);

void Calculate_lnV() {
  const int matSize = ndim*ndim;
  for (int i = 0; i < matSize*nelements; ++i) {
    lnV_avg[i] = 0.0;
  }

  for (int elm = 0; elm < nelements; ++elm) {
    /*The next line says local lnV is part of lnV_avg and starts from
      elm*matSize.So when we assign to lnV, we are in effect assigning
      to lnV_avg.*/
    double *lnV = &lnV_avg[elm*matSize];
    double C_RCauchyGreen[ndim*ndim];
    const int countGP = GaussPoints[elm];
    double preFactor = 1.0/((double)countGP);
    for(int gp = 0; gp < countGP; ++gp) {
      int index = fptr[elm] + matSize * gp;
      double *F_element_gp = &(F[index]);
      // Right Cauchy Greeen: C = F^T * F
      // Right Stretch Tensor: U^2 = C
      // U = sqrt of eigenvalues of C
      // stretchs of U = stretchs of V, the left stretch tensor
      // Abaqus outputs ln(V)
      dgemm_(chy, chn, &ndim, &ndim, &ndim, &preFactor, F_element_gp, &ndim,
             F_element_gp, &ndim, &one, C_RCauchyGreen, &ndim);
    }
    for (int i = 0; i < ndim; i++) {
      //for(int j = 0; j < ndim; j++){
        printf("%f %f %f\n",C_RCauchyGreen[i+0],C_RCauchyGreen[i+1],C_RCauchyGreen[i+2] );
      //}
    }
    printf("\n");

    // Calculate the eigen values of C
    double a, b, c, d, e, f;

    a = C_RCauchyGreen[0]; b = C_RCauchyGreen[1]; c = C_RCauchyGreen[2];
    d = C_RCauchyGreen[4]; e = C_RCauchyGreen[5]; f = C_RCauchyGreen[8];

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

    double eval1, eval2, eval3;
    eval1 = sqrtQ*cos(theta/3.0)+I1;
    eval2 = sqrtQ*cos((theta+2.0*kPi)/3.0)+I1;
    eval3 = sqrtQ*cos((theta+4.0*kPi)/3.0)+I1;

    printf("evals: %f %f %f\n",eval1,eval2,eval3);

  } /* end of loop on elements */


} /* end of function*/
