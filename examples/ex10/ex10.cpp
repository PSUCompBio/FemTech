#include "FemTech.h"
#include "blas.h"

/*Delare Functions*/
void ApplyBoundaryConditions(double Time,double dMax, double tMax);

/* Global Variables/Parameters  - could be moved to parameters.h file?  */
double Time;
int nStep;
int nSteps;
int nPlotSteps = 10;
bool ImplicitStatic = false;
bool ImplicitDynamic = false;
bool ExplicitDynamic = true;
double ExplicitTimeStepReduction = 0.8;
double FailureTimeStep = 1e-11;

int main(int argc, char **argv){

  // Initialize the MPI environment
  MPI_Init(NULL, NULL);
  // Get the number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  // Get the rank of the process
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  if (ReadInputFile(argv[1])) {
    PartitionMesh();
  }

  AllocateArrays();

  if(debug) {
    // Printing local arrays of processor (this section can be removed)
    printf("\neptr array in processor %d after partitioning = ", world_rank);
    for (int i = 0; i <= nelements; i++) {
      printf("%d ", eptr[i]);
    }
    printf("\n");
    printf("\neind array in processor %d after partitioning =", world_rank);
    for (int i = 0; i < nelements; i++) {
      printf(" (%d)  ", i);
      for (int j = eptr[i]; j < eptr[i + 1]; j++) {
        printf("%d ", connectivity[j]);
      }
    }
    printf("\n");
    printf("\nType/PartID of element in processor %d after partitioning = ", world_rank);
    for (int i = 0; i < nelements; i++) {
      printf("%s/%d  ", ElementType[i], pid[i]);
    }
    printf("\n");
    printf("\nSize of coordinates array in processor %d after partitioning = %d\n", world_rank, nnodes * ndim);
    printf("\nCoordinates array in processor %d after partitioning =", world_rank);
    for (int i = 0; i < nnodes; i++) {
      printf(" (%d)  ", i);
      for (int j = 0; j < ndim; j++) {
        printf("%.*f ", 1, coordinates[ndim * i + j]);
      }
    }
    printf("\n");
  }

  /* Write inital, undeformed configuration*/
	Time = 0.0;
	nStep = 0;
  WriteVTU(argv[1], nStep, Time);

  if (ImplicitStatic) {
    // Static solution
    double dMax = 0.1; // max displacment in meters
		double tMax = 1.0;
    ShapeFunctions();
    ReadMaterialProperties();
    ApplyBoundaryConditions(Time,dMax,tMax);
    Assembly((char*)"stiffness");
    ApplySteadyBoundaryConditions();
    SolveSteadyImplicit();
    Time = tMax;
    nStep = 1;
    /* Write final, deformed configuration*/
    WriteVTU(argv[1], nStep, Time);
  }
	else if (ImplicitDynamic) {
    // Dynamic Implicit solution using Newmark's scheme for time integration
    double dt = 0.1;
    double tMax = 10.0;
    double dMax = 0.1;// max displacment in meters
    ShapeFunctions();
    ReadMaterialProperties();
    ApplyBoundaryConditions(Time,dMax,tMax);
    Assembly((char*)"stiffness");
    Assembly((char*)"mass");
    /* beta and gamma of Newmark's scheme */
    double beta = 0.25;
    double gamma = 0.5;
    SolveUnsteadyNewmarkImplicit(beta, gamma, dt, tMax, argv[1]);
  }
	else if (ExplicitDynamic){
		// Dynamic Explcit solution using....
		double dt = 1e-4;
    double tMax = 1e-3; //max simulation time in seconds
    double dMax = 0.1; // max displacment in meters
		double Time=0.0;
		int time_step_counter = 0;
		int plot_counter = 0;

    ShapeFunctions();
    ReadMaterialProperties();
    ApplyBoundaryConditions(Time,dMax,tMax);
    Assembly((char*)"stiffness");
    /*  Step-1: Calculate the mass matrix similar to that of belytschko. */
    Assembly((char*)"mass"); // Add Direct-lumped as an option

	}// end if ExplicitDynamic

	/* Below are things to do at end of program */
	if(world_rank == 0){
		WritePVD(argv[1], nStep, Time);
	}
  FreeArrays();
  MPI_Finalize();
  return 0;
}

void ApplyBoundaryConditions(double Time, double dMax, double tMax){
  double tol = 1e-5;
  int count = 0;
	double AppliedDisp;

	// Apply Ramped Displacment
	if (ExplicitDynamic || ImplicitDynamic){
		AppliedDisp = Time*(dMax/tMax);
  }
  else if (ImplicitStatic){
    AppliedDisp = dMax;
  }


  for(int i=0;i<nnodes;i++){
    // if x value = 0, constrain node to x plane (0-direction)
    if(fabs(coordinates[ndim*i+0]-0.0) <tol){
      boundary[ndim*i+0]=1;
      //printf("node : %d x : %d\n", i, count);
      count = count+1;
    }
    // if y coordinate = 0, constrain node to y plane (1-direction)
    if(fabs(coordinates[ndim*i+1]-0.0) <tol){
      boundary[ndim*i+1]=1;
      //printf("node : %d y : %d\n", i, count);
      count = count+1;
    }
    // if z coordinate = 0, constrain node to z plane (2-direction)
    if(fabs(coordinates[ndim*i+2]-0.0) <tol){
      boundary[ndim*i+2]=1;
      //printf("node : %d z : %d\n", i, count);
      count = count+1;
    }
    // if y coordinate = 1, apply disp. to node = 0.1 (1-direction)
    if(fabs(coordinates[ndim*i+1]-1.0) <tol){
      boundary[ndim*i+1]=1;
      //printf("node : %d y2 : %d\n", i, count);
      count = count+1;
      // note that this may have to be divided into
      // diplacement increments for both implicit and
      // explicit solver. In the future this would be
      // equal to some time dependent function i.e.,
      // CalculateDisplacement to get current increment out
      //  displacment to be applied.
      displacements[ndim*i+1] = AppliedDisp;
    }
  }
  //printf("Time = %3.3e, Applied Disp = %3.3e\n",Time,AppliedDisp);
  return;
}
