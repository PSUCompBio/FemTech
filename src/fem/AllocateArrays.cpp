#include "FemTech.h"

double *displacements;
double *velocities;
double *velocities_half;
double *accelerations;
double *Eavg;
double *fe;
double *fe_prev;
double *fi;
double *fi_prev;
double *f_net;
double *fr_curr;
double *fr_prev;
double *fi_curr;
double *f_damp_prev;
double *f_damp_curr;
double *displacements_prev;
double *accelerations_prev;
double *stepTime;
// Used to compute stress by files in materials folder
double *mat1, *mat2, *mat3, *mat4;

int *boundary;
FILE *energyFile;

void AllocateArrays() {
  const int nDOF = nnodes * ndim;
  // Allocate and initialize global displacements
	displacements=(double*)calloc(nDOF,sizeof(double));
  if (!displacements) {
    printf("ERROR : Error in allocating displacements array\n");
    exit(12);
  }
	//Allocate and initialize global boundary conditions
	boundary=(int*)calloc(nDOF, sizeof(int));
  if (!boundary) {
    printf("ERROR : Error in allocating boundary array\n");
    exit(12);
  }
  if (ImplicitDynamic || ExplicitDynamic) {
    // Velocity and acceleration allocations for unsteady problem
    velocities = (double*)calloc(nDOF, sizeof(double));
    if (!velocities) {
      printf("ERROR : Error in allocating velocities array\n");
      exit(12);
    }
    accelerations = (double*)calloc(nDOF, sizeof(double));
    if (!accelerations) {
      printf("ERROR : Error in allocating accelerations array\n");
      exit(12);
    }
    velocities_half = (double*)calloc(nDOF, sizeof(double));
    if (!velocities_half) {
      printf("ERROR : Error in allocating velocities_half array\n");
      exit(12);
    }
		displacements_prev = (double*)calloc(nDOF, sizeof(double));
    if (!displacements_prev) {
      printf("ERROR : Error in allocating displacements_prev array\n");
      exit(12);
    }
		accelerations_prev = (double*)calloc(nDOF, sizeof(double));
    if (!accelerations_prev) {
      printf("ERROR : Error in allocating accelerations_prev array\n");
      exit(12);
    }
    Eavg = (double*)calloc(nelements*ndim*ndim, sizeof(double));
    if (!Eavg) {
      printf("ERROR : Error in allocating Eavg array\n");
      exit(12);
    }
		fe = (double*)calloc(nDOF, sizeof(double)); // External Nodal force vector
    if (!fe) {
      printf("ERROR : Error in allocating fe array\n");
      exit(12);
    }
		fe_prev = (double*)calloc(nDOF, sizeof(double)); // External Nodal force vector at previous timestep
    if (!fe_prev) {
      printf("ERROR : Error in allocating fe_prev array\n");
      exit(12);
    }
		fi = (double*)calloc(nDOF, sizeof(double)); // Internal Nodal force vector
    if (!fi) {
      printf("ERROR : Error in allocating fi array\n");
      exit(12);
    }
		fi_prev = (double*)calloc(nDOF, sizeof(double)); // Internal Nodal force vector at previous timestep
    if (!fi_prev) {
      printf("ERROR : Error in allocating fi_prev array\n");
      exit(12);
    }
		f_net = (double*)calloc(nDOF, sizeof(double)); // Total Nodal force vector
    if (!f_net) {
      printf("ERROR : Error in allocating f_net array\n");
      exit(12);
    }
		fr_prev = (double*)calloc(nDOF, sizeof(double)); // Reaction Nodal force vector at previous timestep
    if (!fr_prev) {
      printf("ERROR : Error in allocating fr_prev array\n");
      exit(12);
    }
		fr_curr = (double*)calloc(nDOF, sizeof(double));// Reaction Nodal force vector at current timestep
    if (!fr_curr) {
      printf("ERROR : Error in allocating fr_curr array\n");
      exit(12);
    }
		fi_curr = (double*)calloc(nDOF, sizeof(double)); // Internal Nodal force vector at current timestep
    if (!fi_curr) {
      printf("ERROR : Error in allocating fi_curr array\n");
      exit(12);
    }
		f_damp_curr = (double*)calloc(nDOF, sizeof(double)); // Linear Bulk Viscosity Damping Nodal force vector
    if (!f_damp_curr) {
      printf("ERROR : Error in allocating f_damp_curr array\n");
      exit(12);
    }
		f_damp_prev = (double*)calloc(nDOF, sizeof(double)); // Linear Bulk Viscosity Damping Nodal force vector at previous timestep
    if (!f_damp_prev) {
      printf("ERROR : Error in allocating f_damp_prev array\n");
      exit(12);
    }
    energyFile = fopen("energy.dat", "w");
    fprintf(energyFile, "# Energy for FEM\n");
    fprintf(energyFile, "# Time  Winternal   Wexternal   WKE   total\n");

    stepTime = (double*)malloc((MAXPLOTSTEPS)*sizeof(double));
    if (!stepTime) {
      printf("ERROR : Error in allocating stepTime array\n");
      exit(12);
    }
    // Allocations for material temporary computation
    // By default initialize mat1 and mat2
    // If HGO or Linear Elastic material model present, allocate mat3 and mat4
    int matCount = 2;
    for (int i = 0; i < nPIDglobal; ++i) {
      if (materialID[i] == 3) {
        matCount = 4;
      } else {
        if (materialID[i] == 4) {
          matCount = 4;
        }
      }
    }
    const int matSize = ndim*ndim;
    mat1 = (double*)malloc(matSize*sizeof(double));
    mat2 = (double*)malloc(matSize*sizeof(double));
    if (matCount == 4) {
      mat3 = (double*)malloc(matSize*sizeof(double));
      mat4 = (double*)malloc(matSize*sizeof(double));
    }
  }
	return;
}
