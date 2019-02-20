#include "FemTech.h"

double *C;
double rho;

void ReadMaterialProperties() {
	// Create the c matrix with matrix properties
	// TODO(Anil) replace the hard coded material properties with that from input
	// file
	const double E = 1000.0; // TODO(Anil) convert to SI units MPa to Pa?
	const double nu = 0.25;
	rho = 1000.0;
	const double G = E/(2.0*(1.0+nu));
	if (nu == 0.5) {
		printf("ERROR : Incompressible solids does not work with current formulation\n");
	}
	const double c11 = E*(1.0-nu)/((1.0-2.0*nu)*(1.0+nu));
	const double c12 = E*nu/((1.0-2.0*nu)*(1.0+nu));

	// TODO(Anil) relax the assumption of 3D in C Matrix
	// Create C matrix size : 6x6 in colum major format (for blas routines)
	C = (double*)calloc(36, sizeof(double));
	C[0] = C[7] = C[14] = c11;
	C[6] = C[12] = C[13] = c12;
	C[1] = C[2] = C[8] = c12;
	C[21] = C[28] = C[35] = G;
}

int materials() {
   // printf() displays the string inside quotation
   printf("Hello, World and Digital Brain!!\n");
   return 0;
}
