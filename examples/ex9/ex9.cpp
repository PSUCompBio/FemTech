#include "FemTech.h"
#include "blas.h"

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
	//test
  int debug = 0;

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
  WriteVTU(argv[1]);

  // Steady solution
  ShapeFunctions();
  ReadMaterialProperties();
  ReadBoundaryCondition();
  AssembleStiffnessMatrix();
  ApplySteadyBoundaryConditions();
  SolveSteadyImplicit();

  // AssembleMassMatrix();

  FreeArrays();
  MPI_Finalize();
  return 0;
}
