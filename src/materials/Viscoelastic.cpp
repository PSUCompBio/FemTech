#include "FemTech.h"
#include "blas.h"

#include <math.h>
// Viscoelastic material
// Evaluates the PK2 stress tensor

void Viscoelastic(int e, int gp) {
  if (ndim == 2) {
    // 6 values saved per gauss point for 3d
    for (int i = 0; i < 3; i++) {
      int index = pk2ptr[e] + 3 * gp + i;
    }
  }
  if (ndim == 3) {
    const int index = fptr[e] + ndim * ndim * gp;
    const int pide = pid[e];
    const double K = properties[MAXMATPARAMS * pide + 1];    
    const double G0 = properties[MAXMATPARAMS * pide + 2];
    const double Ginf = properties[MAXMATPARAMS * pide + 3];
    const double beta = properties[MAXMATPARAMS * pide + 4];

    const double mu = Ginf + (G0-Ginf)*exp(-beta*Time);
    const double lambda = K - 2.0*mu/3.0;

    // Compute Green-Lagrange Tensor: E= (1/2)*(F^T*F - I)
    double matSize = ndim * ndim;
    double *E = mat1;
    double *F_element_gp = &(F[index]);
    double half = 0.5;
    dgemm_(chy, chn, &ndim, &ndim, &ndim, &half, F_element_gp, &ndim,
           F_element_gp, &ndim, &zero, E, &ndim);
    E[0] -= half;
    E[4] -= half;
    E[8] -= half;

    // Compute 2nd Piola-Kirchhoff Stress
    // S = lambda*tr(E)*I+2*mu*E
    double traceE = E[0] + E[4] + E[8];
    double *S = mat2;
    for (int i = 0; i < matSize; ++i) {
      S[i] = 2.0 * mu * E[i];
    }
    S[0] += lambda * traceE;
    S[4] += lambda * traceE;
    S[8] += lambda * traceE;

    // 6 values saved per gauss point for 3d
    // in voigt notation, sigma11
    pk2[pk2ptr[e] + 6 * gp + 0] = S[0];
    // in voigt notation, sigma22
    pk2[pk2ptr[e] + 6 * gp + 1] = S[4];
    // in voigt notation, sigma33
    pk2[pk2ptr[e] + 6 * gp + 2] = S[8];
    // in voigt notation, sigma23
    pk2[pk2ptr[e] + 6 * gp + 3] = S[7];
    // in voigt notation, sigma13
    pk2[pk2ptr[e] + 6 * gp + 4] = S[6];
    // in voigt notation, sigma12
    pk2[pk2ptr[e] + 6 * gp + 5] = S[3];

    FILE_LOGArraySingle(DEBUGLOGIGNORE, &(pk2[pk2ptr[e]+6*gp]), 6, \
        "Element : %d, GP : %d Stress values", e, gp);
    FILE_LOGMatrix_SINGLE(DEBUGLOGIGNORE, F_element_gp, ndim, ndim, "Printing F Matrix");
    FILE_LOGMatrix_SINGLE(DEBUGLOGIGNORE, E, ndim, ndim, "Printing E Matrix");
    FILE_LOGMatrix_SINGLE(DEBUGLOGIGNORE, S, ndim, ndim, "Printing S Matrix");
  } // if ndim == 3
  return;
}
