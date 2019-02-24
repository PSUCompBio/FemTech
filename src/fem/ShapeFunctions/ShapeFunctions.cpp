#include "FemTech.h"
#include "blas.h"

/* Global Variables */
int *gptr;
// TODO(anil) dsptr[i] = ndim*gptr[i] always. Could be eliminated.
int *dsptr;
int *GaussPoints;
double *shp;
double *dshp;
int *nShapeFunctions;
double *detJacobian;

void ShapeFunctions() {
  // set the debug flag for this file
  int debug = 1;
  // Variables for blas
  char* chy = (char*)"T";
  char* chn = (char*)"N";
  double one = 1.0;
  double zero = 0.0;

  // Create the c matrix with matrix properties
  // TODO(Anil) replace the hard coded material properties with that from input
  // file
  const double E = 1000.0; // TODO(Anil) convert to SI units MPa to Pa?
  const double nu = 0.25;
  const double rho = 1000.0;
  const double G = E/(2.0*(1.0+nu));
  if (nu == 0.5) {
    printf("ERROR : Incompressible solids does not work with current formulation\n");
  }
  const double c11 = E*(1.0-nu)/((1.0-2.0*nu)*(1.0+nu));
  const double c12 = E*nu/((1.0-2.0*nu)*(1.0+nu));
  // TODO(Anil) relax the assumption of 3D in C Matrix
  // Create C matrix size : 6x6 in colum major format (for blas routines)
  double *C = (double*)calloc(36, sizeof(double));
  C[0] = C[7] = C[14] = c11;
  C[6] = C[12] = C[13] = c12;
  C[1] = C[2] = C[8] = c12;
  C[21] = C[28] = C[35] = G;
  int Csize = 6;

  // Global Array - keeps track of how many gauss points there are
  // per element.
  GaussPoints = (int *)malloc(nelements * sizeof(int));
  // Global Array - keeps track of how many shp functions there are
  // per element.
  nShapeFunctions = (int *)malloc(nelements * sizeof(int));

  // Global array -
  // When number of shape functions in an element equals
  // the number of nodes in an element, gptr array looks
  // similar to eptr array. However, the number of quadrature
  // points (a.k.a gauss points) can be different. For example,
  // for a 8-noded hex, there can be 8 gauss points or there
  // could be 1. The gptr array works like the eptr array, but
  // allows this difference to occur.
  gptr = (int *)malloc((nelements+1) * sizeof(int));
  dsptr = (int *)malloc((nelements + 1) * sizeof(int));

  int counter = 0;
  int dshp_counter = 0;

  gptr[0] = 0;
  dsptr[0] = 0;
  for (int i = 0; i < nelements; i++) {
    if (strcmp(ElementType[i], "C3D8") == 0) {
      //GuassPoints per element
      GaussPoints[i] = 8;
      nShapeFunctions[i] = 8;

      // shp function array needs to hold 8
      // shp functions for each of these 8 gauss points
      // for this one element there are 8 gauss points,
      // which each have 8 components in shp function array
      // so for this element I need 8 * 8 positions to hold

      // counter = counter + nShapeFunctions[i]*GaussPoints[i];
      counter += 64;

      // the next counter is used to determine the size of the
      // derivative of shp function, dshp. We expand the slots
      // to account for ndim, since derivatives are taken with
      // respect to chi, eta, and iota.
      // dshp_counter = dshp_counter + (ndim * nShapeFunctions[i]*GaussPoints[i]);
      dshp_counter += 192;
    }
    if (strcmp(ElementType[i], "C3D4") == 0) {
      GaussPoints[i] = 1;
      nShapeFunctions[i] = 4;
      // same argument as above
      // counter = counter + nShapeFunctions[i]*GaussPoints[i];
      // dshp_counter = dshp_counter + (ndim * nShapeFunctions[i]*GaussPoints[i]);
      counter += 4;
      dshp_counter += 12;
    }
    gptr[i + 1] = counter;
    dsptr[i + 1] = dshp_counter;
  }

  // for debugging purposes
  if (debug && 1==1) {
    for (int i = 0; i < nelements; i++) {
      printf("(e.%d) - eptr:[%d->%d] - gptr:[%d->%d] -  dsptr:[%d->%d]\n", i, eptr[i], eptr[i + 1], gptr[i], gptr[i + 1], dsptr[i], dsptr[i + 1]);
    }
    printf("size of shp array = %d \n", counter);
    printf("size of dshp array = %d \n", dshp_counter);
  }

  /*set size of shp array  - this holds shp functions for all elements */
  shp =  (double *)calloc(counter, sizeof(double));
  /*set size of dshp array  - this holds derivatives of shp functions for all elements */
  dshp = (double *)calloc(dshp_counter, sizeof(double));

  /* for debugging */
  if (debug) {
    for (int i = 0; i < nelements; i++) {
      printf("e.%d: int. pts = %d, # shp functions = %d\n", i, GaussPoints[i], nShapeFunctions[i]);
      for (int j = 0; j < GaussPoints[i]; j++) {
        for (int k = 0; k < nShapeFunctions[i]; k++) {
          printf(" %d", gptr[i] + j * GaussPoints[i] + k);
        }
        printf("\n");
      }
      printf("\n");
      for (int j = 0; j < GaussPoints[i]; j++) {
        for (int k = 0; k < nShapeFunctions[i]; k++) {
          for (int l = 0; l < ndim; l++) {
            printf("%d ", dsptr[i]  +  j*GaussPoints[i]*ndim  +  k*ndim + l);
          }
          printf("\n");
        }
        printf("\n");
      }
    }
  }

  for (int i = 0; i < nelements; i++) {
    // Depending on element type call correct shape function library
    // 3D 8-noded hex shape function routine
    if (strcmp(ElementType[i], "C3D8") == 0) {
      double *Chi = (double*)malloc(GaussPoints[i] * ndim * sizeof(double));
      double *GaussWeights = (double*)malloc(GaussPoints[i] * sizeof(double));
      double *detJacobian = (double*)malloc(GaussPoints[i] * sizeof(double));
      GaussQuadrature3D(i, GaussPoints[i], Chi, GaussWeights);
      for (int k = 0; k < GaussPoints[i]; k++) {
        ShapeFunction_C3D8(i, k, Chi, detJacobian);
      }
      // Compute local stiffness
      Stiffness3D(i,detJacobian,GaussWeights);

      // Compute Local Mass matrix
      Mass3D(i,detJacobian,GaussWeights);

      // Free all allocated memories
      free(Chi);
      free(GaussWeights);
      free(detJacobian);
    }

    // 3D 4-noded tet shape function routine
    if (strcmp(ElementType[i], "C3D4") == 0) {
      double *Chi = (double*)malloc(GaussPoints[i]* ndim * sizeof(double));
      double *GaussWeights = (double*)malloc(GaussPoints[i] * sizeof(double));
      GaussQuadrature3D(i, GaussPoints[i], Chi, GaussWeights);
      //printf("chi, eta, iota = %f, %f, %f\n", Chi[0], Chi[1], Chi[2]);
      for (int k = 0; k < GaussPoints[i]; k++) {
        ShapeFunction_C3D4(i, k, GaussPoints[i], Chi);
      }
      free(Chi);
      free(GaussWeights);
    }
  }// loop on nelements

  // for debugging
  if (debug) {
    for (int i = 0; i < nelements; i++) {
      printf("shp array e.%d with %d Gauss points, each with %d shp functions \n", i, GaussPoints[i], nShapeFunctions[i]);
      for (int j = 0; j < GaussPoints[i]; j++) {
        printf("int.%d:\n", j);
        for (int k = 0; k < nShapeFunctions[i]; k++) {
          //printf("%8.5f ", shp[gptr[i] + j * GaussPoints[i] + k]);
          printf(" shp: %4.4f dshp: %8.4f %8.4f %8.4f\n",
              shp[gptr[i] + j * GaussPoints[i] + k],
              dshp[dsptr[i] + j * GaussPoints[i] * ndim + k * ndim + 0],
              dshp[dsptr[i] + j * GaussPoints[i] * ndim + k * ndim + 1],
              dshp[dsptr[i] + j * GaussPoints[i] * ndim + k * ndim + 2]);
        }
        printf("\n");
      }
    }
  }
  return;
}
