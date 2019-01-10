#include "digitalbrain.h"

void updateConnectivityGlobalToLocal(void);
bool PartitionMesh(){

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
  idx_t numflag = 0; // we are using C-style arrays
  idx_t ncommonnodes = 2; // number of nodes elements must have in common
  if(nallelements == 2){
		printf("PartitionMesh.cpp: nallelements = 2 so ncommonnodes set to 8 (for parmetis)\n");
		ncommonnodes = 8;
	}
  const int tpwgts_size = ncon * nparts;
  real_t *tpwgts = (real_t *)malloc(tpwgts_size * sizeof(real_t));
  for (int i = 0; i < tpwgts_size; i++) {
    tpwgts[i] = (real_t)1.0 / (real_t)nparts;
  }

  real_t *ubvec = (real_t *)malloc(ncon * sizeof(real_t));
  for (int i = 0; i < ncon; i++) {
    ubvec[i] = (real_t)1.05; // weight tolerance (same weight tolerance for every weight there is)
  }

  // "elmdist" array describes initial distribution of elements to processors
  // Its size is p + 1, where p is number of processors
  idx_t *elmdist = (idx_t *)calloc(world_size + 1, sizeof(idx_t));
  for (int i = 0, n = 0; i <= world_size; i++) {
    elmdist[i] = n;
    n += (nallelements / world_size);
  }
  elmdist[world_size] += (nallelements % world_size);

  // part array is output of ParMETIS_V3_PartMeshKway() function
  // Size of this array equals to number/count of mesh elements local to processor
  // This array indicates which processor i-th element belongs to
  idx_t edgecut; // output of ParMETIS_V3_PartMeshKway function
  idx_t *part = (idx_t *)calloc(nelements, sizeof(idx_t));

  MPI_Comm comm;
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);
  const int Result = ParMETIS_V3_PartMeshKway(elmdist, eptr, connectivity, elmwgt, \
      &wgtflag, &numflag, &ncon, &ncommonnodes, &nparts, tpwgts, ubvec, options, \
      &edgecut, part, &comm);

  // Output of ParMETIS_V3_PartMeshKway call will be used here
  if (Result == METIS_OK) {
    //printf("\nDebug (proc %d): partArrayFromMETIS\n", world_rank);
    //for(int i = 0; i < nelements; ++i) {
    //  printf("%d\t", part[i]);
    //}
    printf("\n");
    // Using inplace all gatherv to create a local copy of global metis allocation
    int intTask = nallelements/world_size;
    const int from = world_rank*intTask;
    const int to = from + nelements;

    int *partGlobal = (int*)malloc(nallelements*sizeof(int));
    int *elementCountGlobal = (int*)malloc(nallelements*sizeof(int));
    for (int i = 0; i < nelements; ++i) {
      partGlobal[i+from] = (int)part[i];
      elementCountGlobal[i+from] = eptr[i+1]-eptr[i];
    }
    // Get global part array in all processors
    int *counts = (int*)malloc(world_size*sizeof(int));
    for (int i=0; i < world_size; i++) {
      counts[i] = elmdist[i+1]-elmdist[i];
    }
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, partGlobal, counts, elmdist, \
        MPI_INT, MPI_COMM_WORLD);
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, elementCountGlobal, counts, elmdist, \
        MPI_INT, MPI_COMM_WORLD);
    //printf("\nDebug (proc %d): partArray and element count Global\n", world_rank);
    //for(int i = 0; i < nallelements; ++i) {
     // printf("%d\t", partGlobal[i]);
    //}
    //printf("\n");
    //for(int i = 0; i < nallelements; ++i) {
    //  printf("%d\t", elementCountGlobal[i]);
    //}
    //printf("\n");
    // Find new number of elements and number of co-ordinates(not unique)
    int nelementsPart = 0;
    int connectivitySizePart = 0;
    for(int i = 0; i < nallelements; ++i) {
      if (partGlobal[i] == world_rank) {
        nelementsPart = nelementsPart + 1;
        connectivitySizePart += elementCountGlobal[i];
      }
    }
    // Redifining variables to redistribute and store after METIS partitioning
    double *coordinatesPart = (double*)malloc(ndim*connectivitySizePart*sizeof(double));
    int *connectivityPart = (int*)malloc(connectivitySizePart*sizeof(int));
    int *eptrPart = (int*)malloc((nelementsPart+1)*sizeof(int));
    int *pidPart = (int*)malloc(nelementsPart*sizeof(int));
    char **ElementTypePart = (char**)malloc(nelementsPart*sizeof(char*));
    for (int i = 0; i < nelementsPart; ++i) {
      ElementTypePart[i] = (char*)malloc(MAX_ELEMENT_TYPE_SIZE*sizeof(char));
    }
    // Global Part is divided into 3 part 0-from, from-to and 2-nallelements
    // In part 1 receive, in part 2 send/copy if self, in part 3 receive receive
    // Part 2 is the local parts before partitioning
    int currentElementCount = 0;
    eptrPart[0] = 0;
    for (int i = 0; i < from; ++i) {
      if (partGlobal[i] == world_rank) {
        int connectSize = elementCountGlobal[i];
        int location = eptrPart[currentElementCount];
        int source = i/intTask;
        if (source > (world_size-1)) {
          source = world_size-1;
        }
        eptrPart[currentElementCount+1] = location+connectSize;
        MPI_Recv(&(connectivityPart[location]), connectSize, MPI_INT, \
            source, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&(coordinatesPart[ndim*location]), ndim*connectSize, MPI_DOUBLE, \
            source, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&(pidPart[currentElementCount]), 1, MPI_INT, \
            source, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(ElementTypePart[currentElementCount], MAX_ELEMENT_TYPE_SIZE, MPI_CHAR, \
            source, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        currentElementCount = currentElementCount+1;
      }
    }
    for (int i = from, j = 0; i < to; ++i, ++j) {
      if (partGlobal[i] == world_rank) {
        int location = eptr[j];
        int connectSize = eptr[j+1]-location;
        int locationDest = eptrPart[currentElementCount];
        eptrPart[currentElementCount+1] = locationDest+connectSize;
        memcpy(&(connectivityPart[locationDest]), &(connectivity[location]), \
            connectSize*sizeof(int));
        memcpy(&(coordinatesPart[ndim*locationDest]), &(coordinates[ndim*location]), \
            ndim*connectSize*sizeof(double));
        pidPart[currentElementCount] = pid[j];
        memcpy(ElementTypePart[currentElementCount], ElementType[j], \
            MAX_ELEMENT_TYPE_SIZE*sizeof(char));
        currentElementCount = currentElementCount+1;
      } else {
        int dest = partGlobal[i];
        int location = eptr[j];
        int connectSize = eptr[j+1]-location;
        MPI_Send(&(connectivity[location]), connectSize, MPI_INT, \
            dest, i, MPI_COMM_WORLD);
        MPI_Send(&(coordinates[ndim*location]), ndim*connectSize, MPI_DOUBLE, \
            dest, i, MPI_COMM_WORLD);
        MPI_Send(&(pid[j]), 1, MPI_INT, dest, i, MPI_COMM_WORLD);
        MPI_Send(ElementType[j], MAX_ELEMENT_TYPE_SIZE, MPI_CHAR, dest, i, MPI_COMM_WORLD);
      }
    }
    for (int i = to; i < nallelements; ++i) {
      if (partGlobal[i] == world_rank) {
        int connectSize = elementCountGlobal[i];
        int location = eptrPart[currentElementCount];
        int source = i/intTask;
        if (source > (world_size-1)) {
          source = world_size-1;
        }
        eptrPart[currentElementCount+1] = location+connectSize;
        MPI_Recv(&(connectivityPart[location]), connectSize, MPI_INT, \
            source, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&(coordinatesPart[ndim*location]), ndim*connectSize, MPI_DOUBLE, \
            source, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&(pidPart[currentElementCount]), 1, MPI_INT, \
            source, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(ElementTypePart[currentElementCount], MAX_ELEMENT_TYPE_SIZE, MPI_CHAR, \
            source, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        currentElementCount = currentElementCount+1;
      }
    }
    // Swap current variables with variables after partition
    for (int i = 0; i < nelements; ++i) {
      free(ElementType[i]);
    }
    free(ElementType);
    free(connectivity);
    connectivity = connectivityPart;
    free(eptr);
    eptr = eptrPart;
    free(coordinates);
    free(pid);
    pid = pidPart;
    coordinates = coordinatesPart;
    nnodes = eptr[nelementsPart];
    nelements = nelementsPart;
    ElementType = ElementTypePart;
    // Free allocated temporary variables
    free(partGlobal);
    free(counts);

    // Reorder local connectivity
    updateConnectivityGlobalToLocal();
  } else {
    printf("\nERROR( proc %d ): ParMETIS returned error code %d\n", world_rank, Result);
  }
  free(ubvec);
  free(tpwgts);
  free(elmdist);
  free(part);
  return Result == METIS_OK;
}
//-------------------------------------------------------------------------------------------
int compare (const void * a, const void * b) {
  return ( *(int*)a - *(int*)b );
}
//-------------------------------------------------------------------------------------------
int unique(int *arr, int n) {
	int *temp;
	int j = 0;
	temp = (int *)calloc(n, sizeof(int));
  for (int i=0; i < n-1; i++) {
    if (arr[i] != arr[i+1]) {
      temp[j++] = arr[i];
    }
  }
  temp[j++] = arr[n-1];
  for (int i = 0; i < j; i++) {
    arr[i] = temp[i];
  }
	free(temp);
  return j;
}
//-------------------------------------------------------------------------------------------
void updateConnectivityGlobalToLocal(void) {
  int totalSize = eptr[nelements];
  int *newConnectivity = (int*)malloc(totalSize*sizeof(int));
  int *sorted = (int*)malloc(totalSize*sizeof(int));
  memcpy(sorted, connectivity, totalSize*sizeof(int));
  qsort(sorted, totalSize, sizeof(int), compare);
  nnodes = unique(sorted, totalSize);
  for(int i = 0; i < totalSize; ++i) {
    int j;
    for (j = 0; j < nnodes; ++j) {
      if (sorted[j] == connectivity[i]) {
        break;
      }
    }
    newConnectivity[i] = j;
  }
  // Reoder co-ordinates
  double *newCoordinates = (double*)malloc(ndim*nnodes*sizeof(double));
  for (int j = 0; j < nnodes; ++j) {
    int i;
    for(i = 0; i < totalSize; ++i) {
      if (sorted[j] == connectivity[i]) {
        break;
      }
    }
    memcpy(newCoordinates+ndim*j, coordinates+ndim*i, ndim*sizeof(double));
  }
  free(connectivity);
  connectivity = newConnectivity;
  free(coordinates);
  coordinates = newCoordinates;
  free(sorted);
}
