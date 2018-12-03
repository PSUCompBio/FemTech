#include "digitalbrain.h"

//-------------------------------------------------------------------------------------------
int main(int argc, char **argv)
{
    MPI_Initialize();
    
    int ExitStatus = EXIT_FAILURE, partArraySize;
    idx_t *elmdist, *MyEptr, *MyEind, *part, edgecut;
    
    const int NParts = (argc < 3 || atoi(argv[2]) == 0) ? world_size : atoi(argv[2]); // Number of partitions
    
    if (ReadInputFile(argv[1], &elmdist, &MyEptr, &MyEind, &partArraySize))
    {
        if (PartitionMesh(elmdist, MyEptr, MyEind, NParts, partArraySize, &edgecut, &part))
        {
            ExitStatus = EXIT_SUCCESS;
            PARTOUTPUT *Parts = SendReceivePartitioningOutput(edgecut, part, partArraySize);
            if (Parts != NULL)
            {
                // I am processor 0
                CreatePartitions(argv[1], NParts, Parts);
                for (int i = 0; i < world_size; i++)
                {
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
