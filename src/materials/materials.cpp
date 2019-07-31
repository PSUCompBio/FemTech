#include "FemTech.h"

double *C;

void CreateLinearElasticityCMatrix() {
	// Create the c matrix with matrix properties
	const double rho = properties[0];
	const double G = properties[1];
	const double lambda = properties[2];
	const double E = G*(3.0*lambda+2.0*G)/(lambda+G);
	const double nu = 0.5*lambda/(lambda+G);
	if (nu == 0.5) {
		printf("ERROR : Incompressible solids does not work with current formulation\n");
	}
	const double c11 = E*(1.0-nu)/((1.0-2.0*nu)*(1.0+nu));
	const double c12 = E*nu/((1.0-2.0*nu)*(1.0+nu));

	// Create C matrix size : 6x6 in colum major format (for blas routines)
	C = (double*)calloc(36, sizeof(double));
	C[0] = C[7] = C[14] = c11;
	C[6] = C[12] = C[13] = c12;
	C[1] = C[2] = C[8] = c12;
	C[21] = C[28] = C[35] = G;
}
