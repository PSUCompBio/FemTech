#include "FemTech.h"

extern "C" {
 extern double ddot_(int *n, double *dx, int *incx, double *dy, int *incy);
}

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

#if 0
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
#endif

  WriteVTU(argv[1]);


  ShapeFunctions();

	double  m[10],n[10];
	int i;
	double result;
	int len = 10, incm = 1, incn = 1;

	// printf("Enter the elements into first vector.\n");
  result = 0;
	for(i=0;i<10;i++) {
    m[i] = i+1;
    n[i] = i*2;
    result += (i+1)*(i*2);
  }
	printf("The actual result is %f\n",result);
	// scanf("%lf",&m[i]);

	// printf("Enter the elements into second vector.\n");
	// for(i=0;i<10;i++)
	// scanf("%lf",&n[i]);

	result = ddot_(&len, m, &incm, n, &incn);
	printf("The result from blas is %f\n",result);

  // Assembly("mass");

  FreeArrays();
  MPI_Finalize();
  return 0;
}
