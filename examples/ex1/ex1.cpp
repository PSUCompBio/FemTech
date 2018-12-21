#include "digitalbrain.h"

int main(int argc, char **argv)
{
    MPI_Initialize();   
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
