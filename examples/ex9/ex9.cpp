#include "FemTech.h"
#include "blas.h"

/*Delare Functions*/
void ApplyBoundaryConditions();

int main(int argc, char **argv){
	int debug = 0;
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
  WriteVTU(argv[1]);

  // Steady solution
  ShapeFunctions();
  ReadMaterialProperties();
  ReadBoundaryCondition();
  ApplyBoundaryConditions();
  AssembleStiffnessMatrix();
  ApplySteadyBoundaryConditions();
  SolveSteadyImplicit();

  // AssembleMassMatrix();

  FreeArrays();
  MPI_Finalize();
  return 0;
}

void ApplyBoundaryConditions(){
//	0: free
//	1: x prescribed
//	2: y prescribed
//	3: x, y prescribed
//	4: z prescribed
//	5: x, z prescribed
//	6: y, z prescribed
//	7: x, y, z prescribed.

	for(int i=0;i<nnodes;i++){
 		// if x value = 0, constrain node to x plane
		if(fabs(coordinates[ndim*i+0]-0.0) <tol){
			boundary[ndim*i+0]=1;
		}
		// if y coordinate = 0, constrain node to y plane
		if(fabs(coordinates[ndim*i+1]-0.0) <tol){
			boundary[ndim*i+1]=2;
		}
		// if z coordinate = 0, constrain node to z plane
		if(fabs(coordinates[ndim*i+2]-0.0) <tol){
			boundary[ndim*i+2]=4;
		}
		// if y coordinate = 1, apply disp. to node = 0.1
		if(fabs(coordinates[ndim*i+2]-1.0) <tol){
			boundary[ndim*i+1]=2;
			 // note that this may have to be divided into
       // diplacement increments for both implicit and
  		 // explicit solver. In the future this would be
			 // equal to some time dependent function i.e.,
			 // CalculateDisplacement to get current increment out
			 //  displacment to be applied.
			displacements[ndim*i+1] = 0.1
		}
	}



	return;
}
