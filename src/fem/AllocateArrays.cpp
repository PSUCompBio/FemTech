#include "FemTech.h"

double *displacements;
int *materialID;
double *properties;
double *velocities;
double *velocities_half;
double *accelerations;
double *Eavg;
double *fe;
double *fi;
double *f_net;
double *fr_prev;
double *fr_curr;
double *fi_prev;
double *fi_curr;
double *f_damp_prev;
double *f_damp_curr;
double *displacements_prev;


int *boundary;

void AllocateArrays(){

  // Allocate and initialize global displacements
	displacements=(double*)calloc(nnodes*ndim,sizeof(double));
	materialID=(int*)calloc(nelements,sizeof(int));
	properties=(double*)calloc(nelements*MAXMATPARAMS,sizeof(double));

	//Allocate and initialize global boundary conditions
	boundary=(int*)calloc(nnodes*ndim, sizeof(int));
  if (ImplicitDynamic || ExplicitDynamic) {
    // Velocity and acceleration allocations for unsteady problem
    velocities = (double*)calloc(nnodes*ndim, sizeof(double));
    accelerations = (double*)calloc(nnodes*ndim, sizeof(double));
    velocities_half = (double*)calloc(nnodes*ndim, sizeof(double));
		displacements_prev = (double*)calloc(nnodes*ndim, sizeof(double));
    Eavg = (double*)calloc(nelements*ndim*ndim, sizeof(double));
		fe = (double*)calloc(nnodes*ndim, sizeof(double)); // External Nodal force vector
		fi = (double*)calloc(nnodes*ndim, sizeof(double)); // Internal Nodal force vector
		f_net = (double*)calloc(nnodes*ndim, sizeof(double)); // Total Nodal force vector
		fr_prev = (double*)calloc(nnodes*ndim, sizeof(double)); // Reaction Nodal force vector at previous timestep
		fr_curr = (double*)calloc(nnodes*ndim, sizeof(double));// Reaction Nodal force vector at current timestep
		fi_prev = (double*)calloc(nnodes*ndim, sizeof(double)); // Internal Nodal force vector at previous timestep
		fi_curr = (double*)calloc(nnodes*ndim, sizeof(double)); // Internal Nodal force vector at current timestep
		f_damp_curr = (double*)calloc(nnodes*ndim, sizeof(double)); // Linear Bulk Viscosity Damping Nodal force vector
		f_damp_prev = (double*)calloc(nnodes*ndim, sizeof(double)); // Linear Bulk Viscosity Damping Nodal force vector at previous timestep
  }
	return;
}
