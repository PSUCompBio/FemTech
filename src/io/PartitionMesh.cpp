#include "digitalbrain.h"

//-------------------------------------------------------------------------------------------
void MPI_Initialize()
{
    //printf("This is a parallel build!\n");
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    printf("Hello world from processor %s, rank %d out of %d processors\n", processor_name, world_rank, world_size);
}
//-------------------------------------------------------------------------------------------
/*
  This function must be called if/when ReadInputFile() returns true.
  part array is output of this function.
*/
bool PartitionMesh(
    idx_t *elmdist, idx_t *MyEptr, idx_t *MyEind, const idx_t NParts, const int partArraySize, idx_t *edgecut, idx_t **part)
{
    //printf("Let's do some stuff on processor ID %d\n", world_rank);  
    
    // Options for ParMETIS
    idx_t options[3];
    options[0] = 1; // timing info printed
    options[1] = 0; // timing info suppressed
    options[2] = 15; //seed for random number generator

    // Number of partitions (one for each process)
    idx_t nparts = NParts; // number of parts equals number of processes

    // Strange weight arrays needed by ParMETIS
    idx_t ncon = 1; // number of balance constraints

    // Prepare remaining arguments for ParMETIS
    idx_t* elmwgt = NULL;
    idx_t wgtflag = 0; // we don't use weights
    idx_t numflag = 0; // we are using C-style arrays
    idx_t ncommonnodes = 8; // number of nodes elements must have in common
    
    const int tpwgts_size = ncon * nparts;
    real_t *tpwgts = (real_t *)malloc(tpwgts_size * sizeof(real_t));
    for (int i = 0; i < tpwgts_size; i++)
    {
        tpwgts[i] = (real_t)1.0 / (real_t)nparts;
    }
    
    real_t *ubvec = (real_t *)malloc(ncon * sizeof(real_t));
    for (int i = 0; i < ncon; i++)
    {
        ubvec[i] = (real_t)1.05; // weight tolerance (same weight tolerance for every weight there is)
    }

    // part array is output of ParMETIS_V3_PartMeshKway() function
    // Size of this array equals to count of mesh elements local to processor
    // This array indicates which processor i-th element belongs to
    *edgecut = 0; // output of ParMETIS_V3_PartMeshKway function
    *part = (idx_t *)calloc(partArraySize, sizeof(idx_t));
    
    MPI_Comm comm;
    MPI_Comm_dup(MPI_COMM_WORLD, &comm);
    const int Result = ParMETIS_V3_PartMeshKway(
            elmdist,
            MyEptr,
            MyEind,
            elmwgt,
            &wgtflag,
            &numflag,
            &ncon,
            &ncommonnodes,
            &nparts,
            tpwgts,
            ubvec,
            options,
            edgecut,
            *part,
            &comm);

    if (Result != METIS_OK)
    {
        free(*part);
        printf("\nERROR( proc %d ): ParMETIS returned error code %d\n", world_rank, Result);
    }

    free(ubvec);
    free(tpwgts);
    free(elmdist);
    free(MyEptr);
    free(MyEind);
    
    return Result == METIS_OK;
}
