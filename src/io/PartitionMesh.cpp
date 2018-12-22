#include "digitalbrain.h"

bool PartitionMesh(){
	
	//
	// Arrays below need to be defined. They used to be sent in to function
	//
	idx_t *elmdist; //Add comment or description here
	idx_t *MyEptr; //Add comment or description here
	idx_t *MyEind; //Add comment or description here
	idx_t NParts; //Add comment or description here
	int partArraySize; //Add comment or description here
	idx_t *edgecut; //Add comment or description here
	idx_t **part; //Add comment or description here
    
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
