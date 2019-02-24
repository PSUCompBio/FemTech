#include "FemTech.h"

double *displacements;
int *boundary;

void AllocateArrays(){

  // Allocate and initialize global displacements
	displacements=(double*)calloc(nnodes*ndim, sizeof(double));

	//Allocate and initialize global boundary conditions
	boundary=(int*)calloc(nnodes*ndim, sizeof(int));

	return;
}
