#include "FemTech.h"
#include "blas.h"

void Mass3D(int e, double *detJ, double *gweights){
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


	double *Me = (double*)calloc(keSize, sizeof(double));
	double *MeGQ = (double*)calloc(keSize, sizeof(double));
	int nColSize = nShapeFunctions[e]*ndim;
	int Nsize = nColSize*ndim;
	int NiSize = ndim*ndim;
	int NindexStart;
	double Ni;
	// Create N Matrix of size 3x24 in column major format.
	// N is computed for each Gauss Quadrature point
	double *N = (double*)calloc(Nsize, sizeof(double));
	for (int k = 0; k < GaussPoints[e]; k++) {
		// Populate N for each Gauss Point
		for (int n = 0; n < nShapeFunctions[e]; ++n) {
			NindexStart = n*NiSize;
			Ni = shp[gptr[e]+k*GaussPoints[e]+n];

			N[NindexStart] = Ni;
			N[NindexStart+4] = Ni;
			N[NindexStart+8] = Ni;
		}
		// Compute N^T N
		dgemm_(chy, chn, &nColSize, &nColSize, &ndim, &one, N, &ndim, \
				N, &ndim, &zero, MeGQ, &nColSize);
		const double preFactor = gweights[k]*detJ[k];
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
		printf("DEBUG : Printing Me (Mass Matrix) for Element %d\n", e);
		for (int j = 0; j < bColSize; ++j) {
			for (int k = 0; k < bColSize; ++k) {
				printf("%.4f\t", Me[j+k*bColSize]);
			}
			printf("\n");
		}
	}

	free(Me);
	free(MeGQ);
	free(N);

	return;

}
