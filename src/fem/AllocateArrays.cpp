#include "FemTech.h"

double *displacements;
int *boundary;

void AllocateArrays(){

  // Allocate and initialize global displacements
	displacements=(double*)calloc(nnodes*ndim, sizeof(double));
	for(int i=0;i<nnodes;i++){
		for(int j=0;j<ndim;j++){
			displacements[i*ndim+j]=0.0;
		}
	}

	//Allocate and initialize global boundary conditions
	boundary=(int*)calloc(nnodes*ndim, sizeof(int));
	for(int i=0;i<nnodes;i++){
		for(int j=0;j<ndim;j++){
			boundary[i*ndim+j]=0;
		}
	}

	return;
}
