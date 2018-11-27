#include "digitalbrain.h"

//-------------------------------------------------------------------------------------------
// Creates eptr array for processor and returns it size.
// But it doesn't create eind array. MyEind pointer just points to part of the global "connectivity" array.
// Returns size of that local part of the "connectivity" array.
// Returns size of part array
static void GetMyArrays(
    const int world_size,
    const int world_rank,
    const int nelements,
    const idx_t *elmdist,
    idx_t **MyEptr,
    size_t *MyEptrSize,
    idx_t **MyEind,
    size_t *MyEindSize,
    size_t *partArraysize)
{
    const size_t MyElemensCount = elmdist[world_rank + 1] - elmdist[world_rank];
    *MyEptrSize = MyElemensCount + 1;
    *partArraysize = MyElemensCount;
    *MyEindSize = 0;
    *MyEptr = NULL;
    *MyEind = NULL;
    
    int *NodesCountPerElement = (int *)malloc(nelements * sizeof(int));
    for (int i = 1; i <= nelements; i++)
    {
        NodesCountPerElement[i - 1] = eptr[i] - eptr[i - 1];
    }
    
    *MyEptr = (idx_t *)malloc(*MyEptrSize * sizeof(idx_t));

    const size_t From = world_rank * nelements / world_size;
    const size_t To = From + MyElemensCount - 1;
    (*MyEptr)[0] = 0;
    for (int i = 0, j = 0, n = 0, ie = 0; i < nelements; i++)
    {
        if (i >= From && i <= To)
        {
            n += NodesCountPerElement[i];
            (*MyEptr)[++j] = n;
            (*MyEindSize) += NodesCountPerElement[i];
        }
        if ((*MyEind) == NULL && i == From)
        {
            *MyEind = &connectivity[ie];
        }
        ie += NodesCountPerElement[i];
    }
    
    free(NodesCountPerElement);   
}
//-------------------------------------------------------------------------------------------
void PartitionMesh(const char* InpuFileName)
{
    
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
    MPI_Comm comm;
    MPI_Comm_dup(MPI_COMM_WORLD, &comm);

    // Print off a hello world message
    printf("Hello world from processor %s, rank %d out of %d processors\n", processor_name, world_rank, world_size);
#else 
    printf("This is NOT a parallel build!\n");
    int world_rank = 0;
    int world_size = 1;
    MPI_Comm comm = MPI_COMM_WORLD;
#endif

    //read inputfile and initalize
    printf("Let's do some stuff on processor ID %d\n", world_rank);
    printf("This is argv[1]: %s\n", InpuFileName);
    ReadInputFile(InpuFileName);
    
    // Options for ParMETIS
    idx_t options[3];
    options[0] = 1; // timing info printed
    options[1] = 0; // timing info suppressed
    options[2] = 15; //seed for random number generator

    // Number of partitions (one for each process)
    idx_t nparts = world_size; // number of parts equals number of processes

    // Strange weight arrays needed by ParMETIS
    idx_t ncon = 1; // number of balance constraints

    // Prepare remaining arguments for ParMETIS
    idx_t* elmwgt = NULL;
    idx_t wgtflag = 0; // we don't use weights
    idx_t edgecut = 0; // output of ParMETIS_V3_PartMeshKway function
    idx_t numflag = 0; // we are using C-style arrays
    idx_t ncommonnodes = 2; // number of nodes elements must have in common
    
    const size_t tpwgts_size = ncon * nparts;
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

    //ELMDIST: THIS ARRAY DESCRIBES HOW THE ELEMENTS OF THE MESH ARE DISTRIBUTED AMONG THE PROCESSORS.
    //         IT IS ANALOGOUS TO THE VTXDIST ARRAY. ITS CONTENTS ARE IDENTICAL FOR EVERY PROCESSOR.
    //         Size of this array equals to p + 1, where p is count of processors
    idx_t *elmdist = (idx_t *)malloc((world_size + 1) * sizeof(idx_t));
    for (int i = 0, n = 0; i <= world_size; i++)
    {
        elmdist[i] = n;
        n += (nelements / world_size);
    }
    elmdist[world_size] += (nelements % world_size);
     
    printf("\nelmdist array = ");
    for (int i = 0; i <= world_size; i++)
    {
        printf("%d  ", elmdist[i]);
    }

    // Creating local eptr array
    // Size of this array equals to n + 1, where n is count of mesh elements local to processor
    // This array indicates ranges in eind array
    // MyEind points to part of "connectivity" array. Size of this part is returned in MyEindSize.
    // MyEindSize equals to value of MyEptr[MyEptrSize-1]
    idx_t *MyEptr, *MyEind;
    size_t MyEptrSize, MyEindSize, partArraySize;
    GetMyArrays(world_size, world_rank, nelements, elmdist, &MyEptr, &MyEptrSize, &MyEind, &MyEindSize, &partArraySize);
    printf("\nSize of eptr array in processor %d = %d", world_rank, MyEptrSize);
    printf("\neptr array in processor %d = ", world_rank);
    for (size_t i = 0; i < MyEptrSize; i++)
    {
        printf("%d  ", MyEptr[i]);
    }

    printf("\nSize of eind array in processor %d = %d", world_rank, MyEindSize);
    printf("\neind array in processor %d = ", world_rank);
    for (size_t i = 0; i < MyEindSize; i++)
    {
        printf("%d  ", MyEind[i]);
    }
    
    // part array is output of ParMETIS_V3_PartMeshKway() function
    // Size of this array is resolved in GetMyArrays() function and equals to count of mesh elements local to processor
    // This array indicates which processor i-th element belongs to
    printf("\nSize of part array in processor %d = %d", world_rank, partArraySize);
    idx_t *part = (idx_t *)calloc(partArraySize, sizeof(idx_t));
    
    printf("\nCalling ParMETIS_V3_PartMeshKway in processor %d", world_rank);
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
            &edgecut,
            part,
            &comm);

    if (Result == METIS_OK)
    {
        printf("\nPartitioning complete without error in processor %d", world_rank);    
        printf("\nedgecut in processor %d = %d", world_rank, edgecut);
        printf("\npart array in processor %d = ", world_rank);
        for (size_t i = 0; i < partArraySize; i++)
        {
            printf("%d  ", part[i]);
        }
        printf("\n");
    }
    else
    {
        printf("\nParMETIS returned error code %d\n", Result);
    }

    free(part);
    free(MyEptr);
    free(elmdist);
    free(ubvec);
    free(tpwgts);

    WriteVTU((char *)InpuFileName, world_size, world_rank);
    FreeArrays();

#if PARALLEL
    // Finalize the MPI environment.
    MPI_Finalize();
#endif
}
