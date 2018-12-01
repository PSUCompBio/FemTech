#include "digitalbrain.h"

//-------------------------------------------------------------------------------------------
int main(int argc, char **argv)
{
    MPI_Initialize();
    
    int ExitStatus = EXIT_FAILURE;
    idx_t *elmdist, *MyEptr, *MyEind, *part, edgecut;
    size_t partArraySize;
    if (ReadInputFile(argv[1], &elmdist, &MyEptr, &MyEind, &partArraySize))
    {
		// After reading inputfile this block of code should work.
		//printf("nelements = %d\n", nelements);
		//printf("nnodes = %d\n", nnodes);
		//for (int i = 0; i < nelements; i++) {
		//	for (int j = MyEptr[i]; j < MyEptr[i + 1]; j++) {
		//		printf("%d", MyEind[j]);
		//	}
		//}

        if (PartitionMesh(elmdist, MyEptr, MyEind, world_size, partArraySize, &edgecut, &part))
        {
            ExitStatus = EXIT_SUCCESS;
			// After partitioning MyEptr and MyEind are redefined.
			// This block of code should work and be updated.
			//printf("nelements = %d\n", nelements);
			//printf("nnodes = %d\n", nnodes);
			//for (int i = 0; i < nelements; i++) {
			//	for (int j = MyEptr[i]; j < MyEptr[i + 1]; j++) {
			//		printf("%d", MyEind[j]);
			//	}
			//}

            PARTOUTPUT *Parts = SendReceivePartitioningOutput(edgecut, part, partArraySize);
            if (Parts != NULL)
            {
                // I am processor 0.
                for (int i = 0; i < world_size; i++)
                {
                    printf("\nedgecut of processor %d = %d", i, Parts[i].edgecut);
                    printf("\npart array of processor %d = ", i);
                    for (size_t j = 0; j < Parts[i].part_size; j++)
                    {
                        printf("%d ", Parts[i].part[j]);
                    }
                    printf("\n");

					//printf("p%d\n", world_rank);

                    free(Parts[i].part);
                }
                free(Parts);
            }
            free(part);
        }
        free(elmdist);
        free(MyEptr);
        free(MyEind);
    }
    
    MPI_Finalize();
	printf("p%d: nelements=%d\n", world_rank,nelements);
    return ExitStatus;
}
