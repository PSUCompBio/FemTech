#include "digitalbrain.h"

//-------------------------------------------------------------------------------------------
int main(int argc, char **argv)
{
    MPI_Initialize();
    
    int ExitStatus = EXIT_FAILURE, partArraySize;
    idx_t *elmdist, *MyEptr, *MyEind, *part, edgecut;
    
    int NParts = (argc < 3 || atoi(argv[2]) == 0) ? world_size : atoi(argv[2]); // Number of partitions
    
    if (ReadInputFile(argv[1], &elmdist, &MyEptr, &MyEind, &partArraySize))
    {
        if (PartitionMesh(elmdist, MyEptr, MyEind, world_size, partArraySize, &edgecut, &part))
        {
            ExitStatus = EXIT_SUCCESS;
            PARTOUTPUT *Parts = SendReceivePartitioningOutput(edgecut, part, partArraySize);
            if (Parts != NULL)
            {
                // I am processor 0
                PARTITION *Partitions;
                NParts = CreatePartitions(argv[1], NParts, Parts, &Partitions);
                bool Success = true;
                for (int i = 0; Success && i < NParts; i++)
                {
                    Success = PrepareVTUData(argv[1], i, NParts, Partitions[i]);
                    if (Success)
                    {
                        WriteVTU(argv[1], i);
                        FreeArrays();
                    }
                }
                if (NParts > 0)
                {
                    for (int i = 0; i < NParts; i++)
                    {
                        free(Partitions[i].Eptr);
                        free(Partitions[i].Elements);
                    }
                    free(Partitions);                
                }
                if (!Success)
                {
                    ExitStatus = EXIT_FAILURE;
                }
            }
        }
    }
    
    MPI_Finalize();
    
    return ExitStatus;
}
