#include "FemTech.h"
#include "blas.h"

// St. Venant-Kirchhoff
// Evaluates the Cauchy stress tensor

void StVenantKirchhoff(int e, int gp) {
#ifdef DEBUG
  if (debug && 1 == 0) {
    printf("shp array e.%d with %d Gauss points, each with %d shp functions \n",
           e, GaussPoints[e], nShapeFunctions[e]);
    // printf("int.%d:\n", j);
    if (1 == 0) {
      for (int k = 0; k < nShapeFunctions[e]; k++) {
        // printf("%8.5f ", shp[gptr[e] + j * GaussPoints[e] + k]);
        printf(" shp: %4.4f dshp: %8.4f %8.4f %8.4f\n",
               shp[gptr[e] + gp * GaussPoints[e] + k],
               dshp[dsptr[e] + gp * GaussPoints[e] * ndim + k * ndim + 0],
               dshp[dsptr[e] + gp * GaussPoints[e] * ndim + k * ndim + 1],
               dshp[dsptr[e] + gp * GaussPoints[e] * ndim + k * ndim + 2]);
      }
    }
  }
#endif //DEBUG

  // double J = 1;
  // printf("element %d, gauss point %d\n",e,gp);

  if (ndim == 2) {
    // 6 values saved per gauss point for 3d
    for (int i = 0; i < 3; i++) {
      int index = cptr[e] + 3 * gp + i;
      // cauchy[index] = 0.0
      // printf("cauchy[%d]\n",index);
    }
  }
  if (ndim == 3) {
    int index = fptr[e] + ndim * ndim * gp;
    int index2 = detFptr[e] + gp;
    double mu =
        properties[MAXMATPARAMS * e + 1]; // in future will be equal to
                                          // component from  properties array
    double lambda =
        properties[MAXMATPARAMS * e + 2]; // in future will be equal to
                                          // component from  properties array

    // Compute Green-Lagrange Tensor: E= (1/2)*(F^T*F - I)
    double matSize = ndim * ndim;
    double *E = (double *)malloc(matSize * sizeof(double));
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
    double *S = (double *)malloc(matSize * sizeof(double));
    for (int i = 0; i < matSize; ++i) {
      S[i] = 2.0 * mu * E[i];
    }
    S[0] += lambda * traceE;
    S[4] += lambda * traceE;
    S[8] += lambda * traceE;

    // 6 values saved per gauss point for 3d
    // in voigt notation, sigma11
    cauchy[cptr[e] + 6 * gp + 0] = S[0];
    // in voigt notation, sigma22
    cauchy[cptr[e] + 6 * gp + 1] = S[4];
    // in voigt notation, sigma33
    cauchy[cptr[e] + 6 * gp + 2] = S[8];
    // in voigt notation, sigma23
    cauchy[cptr[e] + 6 * gp + 3] = S[7];
    // in voigt notation, sigma13
    cauchy[cptr[e] + 6 * gp + 4] = S[6];
    // in voigt notation, sigma12
    cauchy[cptr[e] + 6 * gp + 5] = S[3];
#ifdef DEBUG
		if(debug){
	    printf("---- Calculate Stress ---\n");
	    printf("Element : %d, GP : %d\n", e, gp);
	    for (int k = 0; k < 6; ++k) {
	      int index_1 = cptr[e]+6*gp+k;
	      printf("Stress value : %12.4E\n", cauchy[index_1]);
	    }
	    printf("---- ----\n");
		}
#endif //DEBUG

#ifdef DEBUG
    if (debug && 1 == 0) {
      printf("Printing F Matrix\n");
      for (int i = 0; i < ndim; ++i) {
        for (int j = 0; j < ndim; ++j) {
          printf("%12.6f\t", F_element_gp[j*ndim+i]);
        }
        printf("\n");
      }
      printf("\n");
      printf("Printing E Matrix\n");
      for (int i = 0; i < ndim; ++i) {
        for (int j = 0; j < ndim; ++j) {
          printf("%12.6f\t", E[j*ndim+i]);
        }
        printf("\n");
      }
      printf("\n");
      printf("Printing S Matrix\n");
      for (int i = 0; i < ndim; ++i) {
        for (int j = 0; j < ndim; ++j) {
          printf("%12.6f\t", S[j*ndim+i]);
        }
        printf("\n");
      }
      printf("\n");
      for (int i = 0; i < 6; i++) {
        int index = cptr[e] + 6 * gp + i;
        printf("cauchy[%d] = %3.3e\n", index, cauchy[index]);
      }
    }
#endif //DEBUG
    free(E);
    free(S);
  } // if ndim == 3
  return;
}
