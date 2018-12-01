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
        if (PartitionMesh(elmdist, MyEptr, MyEind, world_size, partArraySize, &edgecut, &part))
        {
            ExitStatus = EXIT_SUCCESS;
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
    
    return ExitStatus;
}
