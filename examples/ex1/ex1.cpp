#include "digitalbrain.h"

int main(int argc, char **argv){
  
	// Initialize the MPI environment
	MPI_Init(NULL, NULL);

	// Get the number of processes
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	// Get the rank of the process
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	ReadInputFile(argv[1]);
	//PartitionMesh();

	// Printing coordinates
	//printf("\nSize of coordinates array in partition %d = %d\n", i, nnodes * ndim);
	//printf("\nCoordinates array in partition %d =", i);
	//for (int j = 0; j < nnodes; j++){
	//	printf(" (%d)  ", j);
	//	for (int k = 0; k < ndim; k++){
	//		printf("%.*f ", 1, coordinates[ndim * j + k]);
	//	}
	//}
	//printf("\n");

    //WriteVTU(argv[1], i);                    
    FreeArrays();
    MPI_Finalize();
    return 0;
}
