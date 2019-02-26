#include "FemTech.h"
#include "blas.h"

/*Delare Functions*/
void ApplyBoundaryConditions();

int main(int argc, char **argv){
  int debug = 0;
	double TotalTime = 1.0;
	double time_n;
	int step_n;

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
	time_n=0.0;
	step_n=0;
  WriteVTU(argv[1],step_n,time_n);

  // Steady solution
  ShapeFunctions();
  ReadMaterialProperties();
  ApplyBoundaryConditions();
  Assembly((char*)"stiffness");
  ApplySteadyBoundaryConditions();
  SolveSteadyImplicit();

  // Assembly((char*)"mass");

  /* Write final, deformed configuration*/
  time_n=1.0;
  step_n=1;
  WriteVTU(argv[1],step_n,time_n);


	/* Below are things to do at end of program */
	if(world_rank == 0){
		WritePVD(argv[1],step_n,time_n);
	}
  FreeArrays();
  MPI_Finalize();
  return 0;
}

void ApplyBoundaryConditions(){

  double tol = 1e-5;
  int count = 0;
  printf("DEBUG : \n");

  for(int i=0;i<nnodes;i++){
    // if x value = 0, constrain node to x plane (0-direction)
    if(fabs(coordinates[ndim*i+0]-0.0) <tol){
      boundary[ndim*i+0]=1;
      printf("node : %d x : %d\n", i, count);
      count = count+1;
    }
    // if y coordinate = 0, constrain node to y plane (1-direction)
    if(fabs(coordinates[ndim*i+1]-0.0) <tol){
      boundary[ndim*i+1]=1;
      printf("node : %d y : %d\n", i, count);
      count = count+1;
    }
    // if z coordinate = 0, constrain node to z plane (2-direction)
    if(fabs(coordinates[ndim*i+2]-0.0) <tol){
      boundary[ndim*i+2]=1;
      printf("node : %d z : %d\n", i, count);
      count = count+1;
    }
    // if y coordinate = 1, apply disp. to node = 0.1 (1-direction)
    if(fabs(coordinates[ndim*i+1]-1.0) <tol){
      boundary[ndim*i+1]=1;
      printf("node : %d y2 : %d\n", i, count);
      count = count+1;
      // note that this may have to be divided into
      // diplacement increments for both implicit and
      // explicit solver. In the future this would be
      // equal to some time dependent function i.e.,
      // CalculateDisplacement to get current increment out
      //  displacment to be applied.
      displacements[ndim*i+1] = 0.1;
    }
  }
  printf("DEBUG : \n");
  return;
}
