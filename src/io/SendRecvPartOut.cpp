#include "digitalbrain.h"

PARTOUTPUT *SendReceivePartitioningOutput(const idx_t edgecut, idx_t *part, const size_t partArraySize)
{
    PARTOUTPUT *Result = NULL;
    if (world_rank == 0)
    {
        Result = (PARTOUTPUT *)calloc(world_size, sizeof(PARTOUTPUT));
        Result[0].edgecut = edgecut;
        Result[0].part_size = partArraySize;
        Result[0].part = (idx_t *)malloc(partArraySize * sizeof(idx_t));
        memcpy(Result[0].part, part, partArraySize * sizeof(idx_t));
        for (int i = 1; i < world_size; i++)
        {
            memset(&Result[i], 0, sizeof(PARTOUTPUT));
            MPI_Recv(&Result[i].edgecut, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&Result[i].part_size, 1, MPI_UNSIGNED, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            Result[i].part = (idx_t *)malloc(Result[i].part_size * sizeof(idx_t));
            MPI_Recv(Result[i].part, Result[i].part_size, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    else
    {
        MPI_Send(&edgecut, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&partArraySize, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD);
        MPI_Send(part, partArraySize, MPI_INT, 0, 0, MPI_COMM_WORLD);        
    }
    return Result;
}
