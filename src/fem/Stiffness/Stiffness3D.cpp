#include "FemTech.h"
#include "blas.h"

void Stiffness3D(int e, double *detJ, double *gweights){
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

	// for debugging
	if (debug) {
		for (int i = 0; i < nelements; i++) {
			printf("shp array e.%d with %d shp functions\n", i, GaussPoints[i]);
			for (int k = 0; k < nShapeFunctions[i]; k++) {
				printf("int %d | dshp: %8.5f %8.5f %8.5f shp: %8.5f\n", k, dshp[dsptr[i] + k * ndim + 0],
														   dshp[dsptr[i] + k * ndim + 1],
														   dshp[dsptr[i] + k * ndim + 2],
														   shp[gptr[i] + k]);

			}
			printf("\n");
		}
	}

  // number of shape functions * ndim
	int bColSize = nShapeFunctions[e]*ndim;
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
	double *detJacobian = (double*)malloc(GaussPoints[e] * sizeof(double));

	// TODO(Anil) loops can be combined and need to store complete shape functions and
	// derivatives can be eliminated.
	// Computing Ke matrix
	int BiSize = ndim*Csize;
	double dNdx, dNdy, dNdz;
	int indexStart, BindexStart;
	for (int k = 0; k < GaussPoints[e]; k++) {
		// Populate B for each Gauss Point
		for (int n = 0; n < nShapeFunctions[e]; ++n) {
			indexStart = dsptr[e]+(k*GaussPoints[e]+n)*ndim;
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
		const double preFactor = gweights[k]*detJ[k];
		// Ke = \Sum_j w_j (B^T C B Det(J))_j
		for (int n = 0; n < keSize; ++n) {
			Ke[n] += KeGQ[n]*preFactor;
		}
	}

	// print Ke Matrix
	if (debug) {
		printf("RK:DEBUG : Printing Ke (Elemental Stiffness Matrix) for Element %d\n", e);
		for (int j = 0; j < bColSize; ++j) {
			for (int k = 0; k < bColSize; ++k) {
				printf("%.4f\t", Ke[j+k*bColSize]);
			}
			printf("\n");
		}
	}

	free(B);
	free(Ke);
	free(BtC);
	free(KeGQ);

	return;
}
