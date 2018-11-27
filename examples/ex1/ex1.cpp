#include "digitalbrain.h"


//-------------------------------------------------------------------------------------------
// An element that contains a set of nodes
struct Element
{
    std::vector<int> Nodes;

    Element()
    {
    }
    
    Element(const int NodeSize):
        Nodes(NodeSize, 0)
    {
    }
};
//-------------------------------------------------------------------------------------------
void GetMyArrays(
    const int world_size,
    const int world_rank,
    const int nelements,
    const idx_t *elmdist,
    std::vector<int> &My_Eptr,
    std::vector<int> &My_Eind,
    size_t &part_size)
{
    My_Eptr.clear();
    My_Eind.clear();
    part_size = 0;
    
    // Packing all elements read by ReadInputFile() function to All_Elements vector
    // This is just to ease retrieving local eptr and eind arrays from the global variables
    // Actually, local eptr and eind arrays could be read from the input file directly
    std::vector<Element> All_Elements;
    for (size_t i = 1, n = 0; i < (nelements + 1); i++)
    {
        All_Elements.push_back(Element(eptr[i] - eptr[i - 1]));
        std::vector<int> &Nodes = All_Elements.back().Nodes;
        for (size_t j = 0; j < Nodes.size(); j++)
        {
            Nodes[j] = connectivity[n + j];
        }
        n += Nodes.size();
    }

    // Retrieving eptr and eind to int vectors to copy later to local eptr and eind
    const size_t My_Elements_Size = elmdist[world_rank + 1] - elmdist[world_rank];
    part_size = My_Elements_Size;
    const size_t From = world_rank * nelements / world_size;
    const size_t To = From + My_Elements_Size - 1;
    My_Eptr.push_back(0);
    for (size_t i = 0, n = 0; i < All_Elements.size(); i++)
    {
        if (i >= From && i <= To)
        {
            n += All_Elements[i].Nodes.size();
            My_Eptr.push_back(n);
            for (size_t j = 0; j < All_Elements[i].Nodes.size(); j++)
            {
                My_Eind.push_back(All_Elements[i].Nodes[j]);
            }
        }
    }
}
//-------------------------------------------------------------------------------------------
int main(int argc, char **argv)
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
	printf("This is argv[1]: %s\n", argv[1]);
	ReadInputFile(argv[1]);
	
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
	real_t *tpwgts = new real_t[tpwgts_size];
	for (int i = 0; i < tpwgts_size; i++)
	{
	    tpwgts[i] = (real_t)1.0 / (real_t)nparts;
	}
	
	real_t *ubvec = new real_t[ncon];
	for (int i = 0; i < ncon; i++)
	{
	    ubvec[i] = (real_t)1.05; // weight tolerance (same weight tolerance for every weight there is)
	}

    //ELMDIST: THIS ARRAY DESCRIBES HOW THE ELEMENTS OF THE MESH ARE DISTRIBUTED AMONG THE PROCESSORS.
	//         IT IS ANALOGOUS TO THE VTXDIST ARRAY. ITS CONTENTS ARE IDENTICAL FOR EVERY PROCESSOR.
	//         Size of this array equals to p + 1, where p is count of processors
 	idx_t *elmdist = new idx_t[world_size + 1];
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

    // Retrieving eptr and eind to int vectors to copy later to local eptr and eind
    std::vector<int> My_Eptr_V, My_Eind_V;
    size_t part_size;
    GetMyArrays(world_size, world_rank, nelements, elmdist, My_Eptr_V, My_Eind_V, part_size);
    
    // Initializing local eptr array
    // Size of this array equals to n + 1, where n is count of mesh elements local to processor
    // This array indicates ranges in eind array
    idx_t * my_eptr = new idx_t[My_Eptr_V.size()];
    printf("\nSize of eptr array in processor %d = %d", world_rank, My_Eptr_V.size());
    printf("\neptr array in processor %d = ", world_rank);
    for (size_t i = 0; i < My_Eptr_V.size(); i++)
    {
        my_eptr[i] = My_Eptr_V[i];
        printf("%d  ", my_eptr[i]);
    }

    // Initializing local eind array
    // Size of this array equals to value of eptr[size_of_eptr-1], where eptr[size_of_eptr-1] is value of the last element of eptr array
    // This array stores nodes of elements that are local to processor
    idx_t * my_eind = new idx_t[My_Eind_V.size()];
    printf("\nSize of eind array in processor %d = %d", world_rank, My_Eind_V.size());
    printf("\neind array in processor %d = ", world_rank);
    for (size_t i = 0; i < My_Eind_V.size(); i++)
    {
        my_eind[i] = My_Eind_V[i];
        printf("%d  ", my_eind[i]);
    }
	
	// part array is output of ParMETIS_V3_PartMeshKway() function
	// Size of this array is resolved in GetMyArrays() function and equals to count of mesh elements local to processor
	// This array indicates which processor i-th element belongs to
    printf("\nSize of part array in processor %d = %d", world_rank, part_size);
	idx_t *part = new idx_t[part_size];
    for (size_t i = 0; i < part_size; i++)
	{
	    part[i] = 0;
	}
	
	//WriteVTU(argv[1]);
	//exit(0);


    //int ParMETIS_V32_Mesh2Dual(idx_t *elmdist, idx_t *eptr, idx_t *eind, idx_t *numflag, idx_t *ncommonnodes, 
		//idx_t **xadj, idx_t **adjncy, MPI_Comm *comm);

    printf("\nCalling ParMETIS_V3_PartMeshKway in processor %d", world_rank);
	const int Result = ParMETIS_V3_PartMeshKway(
	        elmdist,
	        my_eptr,
	        my_eind,
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
		for (size_t i = 0; i < part_size; i++)
		{
		    printf("%d  ", part[i]);
		}
		printf("\n");
	}
	else
	{
		printf("\nParMETIS returned error code %d\n", Result);
	}

    delete [] part;
    delete [] my_eind;    
    delete [] my_eptr;
    delete [] elmdist;
    delete [] ubvec;
    delete [] tpwgts;
    

	//deallocate(options)
	//deallocate(ubvec)
	//deallocate(tpwgts)
	//deallocate(wgt)
	//deallocate(eptr)
	//deallocate(elmdist)
	//deallocate(nelempa)

	WriteVTU(argv[1],world_size, world_rank);
    FreeArrays();

#if PARALLEL
	// Finalize the MPI environment.
	MPI_Finalize();
#endif

    return 0;
}
