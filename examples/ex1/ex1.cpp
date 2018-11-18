#include "digitalbrain.h"

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

	//if (world_rank == 0) {
		//read inputfile and initalize
		printf("Let's do some stuff on processor ID %d\n", world_rank);
		printf("this is argv[1]: %s\n", argv[1]);
		ReadInputFile(argv[1]);
	//}
#if 1
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
	idx_t edgecut = 0;
	idx_t numflag = 0; // we are using C-style arrays
	idx_t ncommonnodes = 2; // number of nodes elements must have in common

	//real_t tpwgts[4];

	real_t *tpwgts = (real_t*)malloc((ncon*nparts) * sizeof(real_t));
	for (int i = 0; i < ncon*nparts; i++) { tpwgts[i] = (real_t)1.0 / (real_t)nparts; }
	
	real_t *ubvec = (real_t*)malloc((ncon) * sizeof(real_t));
	for (int i = 0; i < ncon; i++) { ubvec[i] = (real_t)1.05; } // weight tolerance (same weight tolerance for every weight there is)

    //ELMDIST: THIS ARRAY DESCRIBES HOW THE ELEMENTS OF THE MESH ARE DISTRIBUTED AMONG THE PROCESSORS.
	//         IT IS ANALOGOUS TO THE VTXDIST ARRAY.ITS CONTENTS ARE IDENTICAL FOR EVERY PROCESSOR.
	idx_t *elmdist = (idx_t*)malloc((nparts+1) * sizeof(idx_t));
	elmdist[0] = 0;
	for (int i = 1; i < nparts+1; i++) { 
		//elmdist[i] = elmdist[i-1] + nelements/nparts;
		elmdist[i] = 0; // all elements are on process zero
		printf("elmdist[%d] = %d\n", i, elmdist[i]);
	}

	idx_t *part = (idx_t*)malloc((nelements) * sizeof(idx_t));
	for (int i = 0; i < nelements; i++) { part[i] = 0; }
#endif

	// Create and fill arrays "eptr", where eptr[i] is the number of vertices that belong to the i-th element, and
	// "eind" contains the vertex-numbers of the i-the element in eind[eptr[i]] to eind[eptr[i+1]-1]
	for (int i = 0; i < nelements; i++) {
		for (int j = eptr[i]; j < eptr[i+1]; j++) {
			//printf("%d %d %d \n",i,j, connectivity[j] );

		}
		//printf("\n");
	}
	//WriteVTU(argv[1]);
	//exit(0);

	printf("size of eind = %d\n", eptr[nelements]);
	idx_t *eind= (idx_t*)malloc(eptr[nelements] * sizeof(idx_t));
	for (int i = 0; i < eptr[nelements]; i++) {
		eind[i] = connectivity[i];
	}


	//if (world_rank == 0) {
	    int ParMETIS_V32_Mesh2Dual(idx_t *elmdist, idx_t *eptr, idx_t *eind, idx_t *numflag, idx_t *ncommonnodes, 
			idx_t **xadj, idx_t **adjncy, MPI_Comm *comm);

	

		int ParMETIS_V3_PartMeshKway(idx_t *elmdist, idx_t *eptr, idx_t *eind, idx_t *elmwgt,
			idx_t *wgtflag, idx_t *numflag, idx_t *ncon, idx_t *ncommonnodes, idx_t *nparts,
			real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *edgecut, idx_t *part, MPI_Comm *comm);

		if (METIS_OK) {
			printf("Partitioning complete!\n");
		}
		else {
			printf("ParMETIS returned error code %d\n", METIS_OK);
		}
	//}
	for (int i = 0; i < nelements; i++) {
		printf("part: %d\n", part[i]);
	}
	for (int i = 0; i < nparts+1; i++) {
		printf("elmdist: %d\n", elmdist[i]);
	}

	//deallocate(options)
	//deallocate(ubvec)
	//deallocate(tpwgts)
	//deallocate(wgt)
	//deallocate(eptr)
	//deallocate(elmdist)
	//deallocate(nelempa)




   
   //FreeArrays();

#if PARALLEL
	// Finalize the MPI environment.
	MPI_Finalize();
#endif
   return 0;
}
