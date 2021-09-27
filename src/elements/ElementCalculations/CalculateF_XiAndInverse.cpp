#include "FemTech.h"
#include "blas.h"

// Function to compute F_Xi and its inverse along with the determinant
// The function is used for updated Lagrangian internal force computation
// x_i = x_iI N_I(Xi_i), dx_i/dXi_i = x_iI dN_I/dXi_i
// Assumes dshp stores dN_I/dXi_i

void CalculateF_XiAndInverse(int e, int gp) {
  // if ((eptr[e + 1] - eptr[e]) != nShapeFunctions[e]) {
  //   FILE_LOG_SINGLE(WARNING, "CalculateDeformationGradient.cpp: (eptr[e+1]-eptr[e]) "
  //          "!= nShapeFunctions[e])\nCheck it out, probally a bug in allocation");
  // }
  const unsigned int start = eptr[e];
  const unsigned int nShapeFunc = eptr[e+1]-start;
  const unsigned int indexStart = dsptr[e] + gp * GaussPoints[e] * ndim;
  for (unsigned int i = 0; i < ndim; i++) {
    for (unsigned int j = 0; j < ndim; j++) {
      double F_Xi_ij = 0.0;
      for (unsigned int I = 0; I < nShapeFunc; ++I) {
        const unsigned int dispIndex = ndim*connectivity[start + I] + i;
        const unsigned int indexShp = indexStart + I * ndim + j;
        F_Xi_ij = F_Xi_ij + (coordinates[dispIndex]+displacements[dispIndex])*dshp[indexShp];
      } // loop on I
      // Store in column major format
      F_Xi[i+j*ndim] = F_Xi_ij;
    }
  }

  FILE_LOGMatrix_SINGLE(DEBUGLOGIGNORE, F_Xi, ndim, \
      ndim, "--------F_Xi for gauss point %d --------", gp);
  // double *FTest = (double*)malloc(ndim2*sizeof(double));
  // const int indexT = fptr[e] + ndim * ndim * gp;
  // double *F_element_gp = &(F[indexT]);
  // double *F_Xi_0_egp = &(F_Xi_0[indexT]);
  // dgemm_(chn, chn, &ndim, &ndim, &ndim, &one, F_element_gp, &ndim,
  //          F_Xi_0_egp, &ndim, &zero, FTest, &ndim);
  // for (unsigned int i = 0; i < ndim2; ++i) {
  //   FTest[i] -= F_Xi[i];
  // }
  // FILE_LOGMatrix_SINGLE(WARNING, FTest, ndim, \
  //     ndim, "--------F_Test for gauss point %d --------", gp);
  // FILE_LOGMatrix_SINGLE(WARNING, F_element_gp, ndim, \
  //     ndim, "--------F for gauss point %d --------", gp);
  // FILE_LOGMatrix_SINGLE(WARNING, F_Xi_0_egp, ndim, \
  //     ndim, "--------F_Xi_0 for gauss point %d --------", gp);
  // FILE_LOGMatrix_SINGLE(WARNING, F_Xi, ndim, \
  //     ndim, "--------F_Xi for gauss point %d --------", gp);
  
  // Compute the matrix inverse
  J_Xi = inverse3x3Matrix(F_Xi, F_XiInverse);

  // dgemm_(chn, chn, &ndim, &ndim, &ndim, &one, F_Xi, &ndim,
  //          F_XiInverse, &ndim, &zero, FTest, &ndim);
  // FILE_LOGMatrix_SINGLE(WARNING, FTest, ndim, \
  //     ndim, "--------F_Xi*Inverse for gauss point %d --------", gp);
  // free(FTest);
  return;
}
