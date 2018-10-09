#include "digitalbrain.h"
//#if PARALLEL
//#  define PARALLEL @PARALLEL@
//#endif

int main(int argc, char **argv){
	

#if PARALLEL
		printf("This is a parallel build!\n");
		// Initialize the MPI environment
		MPI_Init(NULL, NULL);

		// Get the number of processes
		int world_size;
		MPI_Comm_size(MPI_COMM_WORLD, &world_size);

		// Get the rank of the process
		int world_rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

		// Get the name of the processor
		char processor_name[MPI_MAX_PROCESSOR_NAME];
		int name_len;
		MPI_Get_processor_name(processor_name, &name_len);

		// Print off a hello world message
		printf("Hello world from processor %s, rank %d out of %d processors\n",
			processor_name, world_rank, world_size);

#else 
		printf("This is NOT a parallel build!\n");
		int world_rank = 0;
#endif

	if (world_rank == 0) {
			//read inputfile and initalize
			printf("Let's do some stuff on processor ID %d\n", world_rank);
			printf("this is argv[1]: %s\n", argv[1]);
			ReadInputFile(argv[1]);
	}

   //WriteVTU(argv[1]);
   //FreeArrays();

#if PARALLEL
	   // Finalize the MPI environment.
	   MPI_Finalize();
   
#endif
   return 0;
}
