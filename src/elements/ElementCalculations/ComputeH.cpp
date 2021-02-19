#include "FemTech.h"
#include "blas.h"

// Function to compute H matrix. H is  used to accurately compute F and E
// matrices in updated Lagrangian Formulation
// H_ij = dN_I/dX_j u_iI

void ComputeH(int e, int gp, double *H) {
  const unsigned int start = eptr[e];
  const unsigned int nShapeFunc = eptr[e+1]-start;
  // const unsigned int nShapeFunc = nShapeFunctions[e];
  const unsigned int indexStart = dsptr[e] + gp * GaussPoints[e] * ndim;
  for (unsigned int i = 0; i < ndim; ++i) {
    for (unsigned int j = 0; j < ndim; ++j) {
      double Hij = 0.0;
      for (unsigned int I = 0; I < nShapeFunc; ++I) {
        const unsigned int dispIndex = connectivity[start+I]*ndim + i;
        const unsigned int indexShp = indexStart + I * ndim + j;
        Hij = Hij + B0[indexShp]*displacements[dispIndex];
      }
      // Store Hij in column major format
      H[i+j*ndim] = Hij;
    }
  }
  // double *FTest = (double*)malloc(ndim2*sizeof(double));
  // const int indexT = fptr[e] + ndim * ndim * gp;
  // double *F_element_gp = &(F[indexT]);
  // for (unsigned int i = 0; i < ndim2; ++i) {
  //   FTest[i] = H[i];
  // }
  // FTest[0] = FTest[0] + 1.0;
  // FTest[4] = FTest[4] + 1.0;
  // FTest[8] = FTest[8] + 1.0;
  // for (unsigned int i = 0; i < ndim2; ++i) {
  //   FTest[i] -= F_element_gp[i];
  // }
  // FILE_LOGMatrix_SINGLE(WARNING, FTest, ndim, \
  //     ndim, "--------H_Test for gauss point %d --------", gp);
}
