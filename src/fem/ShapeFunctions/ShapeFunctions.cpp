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
      const int nShapeFunctions = 8;
      const int nGaussPoints = GaussPoints[i];
      // TODO(Anil) All the below allocations can be hoisted out of loop and 
      // size can be based on the maximum.
      int bColSize = nShapeFunctions*ndim;
      int Bsize = bColSize*Csize;
      int keSize = bColSize*bColSize;
      // Create B Matrix of size 24x6 in column major format.
      // B is computed for each Gauss Quadrature point
      double *B = (double*)calloc(Bsize, sizeof(double));
      // Create local element stiffness matrix Ke of size 24x24 in column major
      // format
      double *Ke = (double*)calloc(keSize, sizeof(double));
      // Create temperory variables to store intermediate results
      double *BtC = (double*)calloc(Bsize, sizeof(double));
      double *KeGQ = (double*)calloc(keSize, sizeof(double));

      double *Chi = (double*)malloc(nGaussPoints * ndim * sizeof(double));
      double *GaussWeights = (double*)malloc(nGaussPoints * sizeof(double));
      double *detJacobian = (double*)malloc(nGaussPoints * sizeof(double));
      GaussQuadrature3D(i, nGaussPoints, Chi, GaussWeights);
      for (int k = 0; k < nGaussPoints; k++) {
        ShapeFunction_C3D8(i, k, Chi, detJacobian);
      }
      // TODO(Anil) loops can be combined and need to store complete shape functions and
      // derivatives can be eliminated.
      // Computing Ke matrix
      int BiSize = ndim*Csize;
      double dNdx, dNdy, dNdz;
      int indexStart, BindexStart;
      for (int k = 0; k < nGaussPoints; k++) {
        // Populate B for each Gauss Point
        for (int n = 0; n < nShapeFunctions; ++n) {
          indexStart = dsptr[i]+(k*nGaussPoints+n)*ndim;
          BindexStart = n*BiSize;
          dNdx = dshp[indexStart];
          dNdy = dshp[indexStart+1];
          dNdz = dshp[indexStart+2];

          B[BindexStart] = dNdx;
          B[BindexStart+7] = dNdy;
          B[BindexStart+14] = dNdz;

          B[BindexStart+9] = dNdz;
          B[BindexStart+15] = dNdy;
          B[BindexStart+4] = dNdz;
          B[BindexStart+16] = dNdx;
          B[BindexStart+5] = dNdy;
          B[BindexStart+11] = dNdx;
        }
        // Compute B^T C
        dgemm_(chy, chn, &bColSize, &Csize, &Csize, &one, B, &Csize, \
            C, &Csize, &zero, BtC, &bColSize);
        // Compute B^T C B
        dgemm_(chn, chn, &bColSize, &bColSize, &Csize, &one, BtC, &bColSize, \
            B, &Csize, &zero, KeGQ, &bColSize);
        const double preFactor = GaussWeights[k]*detJacobian[k];
        // Ke = \Sum_j w_j (B^T C B Det(J))_j
        for (int n = 0; n < keSize; ++n) {
          Ke[n] += KeGQ[n]*preFactor;
        }
      }
      // print Ke Matrix
      if (debug) {
        printf("DEBUG : Printing Ke (Elemental Stiffness Matrix) for Element %d\n", i);
        for (int j = 0; j < bColSize; ++j) {
          for (int k = 0; k < bColSize; ++k) {
            printf("%.4f\t", Ke[j+k*bColSize]);
          }
          printf("\n");
        }
      }
      // Computing Me matrix
      // Create local element mass matrix Me of size 24x24 in column major
      // format
      double *Me = (double*)calloc(keSize, sizeof(double));
      double *MeGQ = (double*)calloc(keSize, sizeof(double));
      int nColSize = nShapeFunctions*ndim;
      int Nsize = nColSize*ndim;
      int NiSize = ndim*ndim;
      int NindexStart;
      double Ni;
      // Create N Matrix of size 3x24 in column major format.
      // N is computed for each Gauss Quadrature point
      double *N = (double*)calloc(Nsize, sizeof(double));
      for (int k = 0; k < nGaussPoints; k++) {
        // Populate N for each Gauss Point
        for (int n = 0; n < nShapeFunctions; ++n) {
          NindexStart = n*NiSize;
          Ni = shp[gptr[i]+k*nGaussPoints+n];

          N[NindexStart] = Ni;
          N[NindexStart+4] = Ni;
          N[NindexStart+8] = Ni;
        }
        // Compute N^T N
        dgemm_(chy, chn, &nColSize, &nColSize, &ndim, &one, N, &ndim, \
            N, &ndim, &zero, MeGQ, &nColSize);
        const double preFactor = GaussWeights[k]*detJacobian[k];
        // Me = \Sum_j w_j (N^T N Det(J))_j
        for (int n = 0; n < keSize; ++n) {
          Me[n] += MeGQ[n]*preFactor;
        }
      }
      printf("\n\nRho : %.4f\n", rho);
      for (int n = 0; n < keSize; ++n) {
        Me[n] *= rho;
      }
      // print Me Matrix
      if (debug) {
        printf("DEBUG : Printing Me (Mass Matrix) for Element %d\n", i);
        for (int j = 0; j < bColSize; ++j) {
          for (int k = 0; k < bColSize; ++k) {
            printf("%.4f\t", Me[j+k*bColSize]);
          }
          printf("\n");
        }
      }
      // Free all allocated memories
      free(Chi);
      free(GaussWeights);
      free(detJacobian);
      free(B);
      free(Ke);
      free(BtC);
      free(KeGQ);
      free(Me);
      free(MeGQ);
      free(N);
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
