#include "FemTech.h"

double *displacements;
double *velocities;
double *accelerations;
int *boundary;

void AllocateArrays(){

  // Allocate and initialize global displacements
	displacements=(double*)calloc(nnodes*ndim, sizeof(double));

	//Allocate and initialize global boundary conditions
	boundary=(int*)calloc(nnodes*ndim, sizeof(int));
  if (unsteadyFlag) {
    // Velocity and acceleration allocations for unsteady problem
    velocities = (double*)calloc(nnodes*ndim, sizeof(double));
    accelerations = (double*)calloc(nnodes*ndim, sizeof(double));
  }
	return;
}
