#include "digitalbrain.h"

//-------------------------------------------------------------------------------------------
int main(int argc, char **argv)
{
    if (argv[1] == NULL)
    {
        printf("Input filename is empty. Quiting.\n");
        return EXIT_FAILURE;
    }
    
    MPI_Initialize();
    
    int ExitStatus = EXIT_SUCCESS;
    if (ReadInputFile(argv[1]))
    {
        PartitionMesh();
        WriteVTU(argv[1]);
        FreeArrays();
    }
    else
    {
        printf("\nReadInputFile() failed.\n");
        ExitStatus = EXIT_FAILURE;
    }
    
    MPI_Finalize();
    
    return ExitStatus;
}
