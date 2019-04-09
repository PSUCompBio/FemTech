#include "FemTech.h"

#include <assert.h>

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
    // Variable to count number of elements to move to each MPI process
    // from the current local set of elements
    // variable part the output of parmetis is used to populate this array
    int *elementRedistributePattern = (int*)calloc(world_size, sizeof(int));
    int *nodeRedistributePattern = (int*)calloc(world_size, sizeof(int));
    int processorToMoveElement, cordInElement;
    for (int i = 0; i < nelements; ++i) {
      processorToMoveElement = (int)part[i];
      cordInElement = eptr[i+1]-eptr[i];
      elementRedistributePattern[processorToMoveElement] += 1;
      nodeRedistributePattern[processorToMoveElement] += cordInElement;
    }
    // Create cumulative array to act as sendbuffer pointer
    int *elementSendCum = (int*)malloc((world_size+1)*sizeof(int));
    int *nodeSendCum = (int*)malloc((world_size+1)*sizeof(int));
    elementSendCum[0] = 0;
    nodeSendCum[0] = 0;
    // Variables to keep track of requests
    int sendCount = 0, recvCount = 0;
    for (int i = 0; i < world_size; ++i) {
      if (elementRedistributePattern[i] && i != world_rank) {
        sendCount += 1;
      }
      elementSendCum[i+1] = elementSendCum[i]+elementRedistributePattern[i];
      nodeSendCum[i+1] = nodeSendCum[i]+nodeRedistributePattern[i];
    }
    // Multiply by 5 as we are sending 5 seperate lists for each process
    sendCount*=5;
    int nnodes_process = nodeSendCum[world_size];
    // Verify that all elements are included for distribution
    assert(nelements == elementSendCum[world_size]);

    // Get count to receive from each processor
    int *elementReceivePattern = (int*)malloc(world_size*sizeof(int));
    int *nodeReceivePattern = (int*)malloc(world_size*sizeof(int));
    memcpy(elementReceivePattern, elementRedistributePattern, world_size*sizeof(int));
    memcpy(nodeReceivePattern, nodeRedistributePattern, world_size*sizeof(int));
    MPI_Alltoall(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, elementReceivePattern, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Alltoall(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, nodeReceivePattern, 1, MPI_INT, MPI_COMM_WORLD);
    if (debug && 1==0) {
      printf("\nDEBUG (proc %d): Element Redistribution Pattern\n", world_rank);
      for (int i = 0; i < world_size; ++i) {
        printf("%d\t",elementRedistributePattern[i]);
      }
      printf("\nDEBUG (proc %d): Element Receive Pattern\n", world_rank);
      for (int i = 0; i < world_size; ++i) {
        printf("%d\t",elementReceivePattern[i]);
      }
      printf("\nDEBUG (proc %d): Node Redistribution Pattern\n", world_rank);
      for (int i = 0; i < world_size; ++i) {
        printf("%d\t",nodeRedistributePattern[i]);
      }
      printf("\nDEBUG (proc %d): Node Receive Pattern\n", world_rank);
      for (int i = 0; i < world_size; ++i) {
        printf("%d\t",nodeReceivePattern[i]);
      }
      printf("\n");
    }
    // Pack data to redistribute mesh
    // Data to send and receive include node coordinates, element connectivity,
    // eptr array, part id array, Element type array
    // Create send buffers

    // Create array to keep track of current count during packing
    int *elementCount = (int*)calloc(world_size, sizeof(int));
    int *nodeCount = (int*)calloc(world_size, sizeof(int));

    // Variables to reorder local data to send
    int *connectivityPacked = (int*)malloc(nnodes_process*sizeof(int));
    double *coordinatesPacked = (double*)malloc(nnodes_process*sizeof(double)*ndim);
    int *eptrPacked = (int*)malloc(nelements*sizeof(int));
    int *pidPacked = (int*)malloc(nelements*sizeof(int));
    char *ElementTypePacked = (char*)malloc(nelements*MAX_ELEMENT_TYPE_SIZE*sizeof(char*));

    // Move data to packed arrays
    for (int i = 0; i < nelements; ++i) {
      // Anil(TODO) This can be optimized to move multiple elements if
      // concecutive elements to move are assigned same destinaation process
      int start = eptr[i];
      int end = eptr[i+1];
      processorToMoveElement = (int)part[i];
      int count = end-start;
      int locationToCopy = nodeSendCum[processorToMoveElement]+\
                           nodeCount[processorToMoveElement]; 
      int elementLocToCopy = elementSendCum[processorToMoveElement]+\
                             elementCount[processorToMoveElement];
      // Copy connectivity 
      memcpy(&(connectivityPacked[locationToCopy]), &(connectivity[start]), \
          sizeof(int)*count);
      // Copy coordinates
      memcpy(&(coordinatesPacked[locationToCopy*ndim]), &(coordinates[start*ndim]), \
          sizeof(double)*count*ndim);
      // Store number of nodes to recreate eptr on destination
      eptrPacked[elementLocToCopy] = count;
      // Copy pid to send
      pidPacked[elementLocToCopy] = pid[i];
      // Copy element type
      memcpy(&(ElementTypePacked[elementLocToCopy*MAX_ELEMENT_TYPE_SIZE]), 
          ElementType[i], sizeof(char)*MAX_ELEMENT_TYPE_SIZE);

      // Advance next write memory location
      nodeCount[processorToMoveElement] += count;
      elementCount[processorToMoveElement] += 1;
    }
    // Verify all elements have been moved
    for (int i = 0; i < world_size; ++i) {
      assert(nodeCount[i] == nodeRedistributePattern[i]);
      assert(elementCount[i] == elementRedistributePattern[i]);
    }

    // Create MPI non-blocking sends of the packed data
    // One MPI send is used per process pair by using the packed structure above
    // Create an array of requests to keep track of completion of the send calls
    MPI_Request* requestListSend = (MPI_Request*)malloc(sizeof(MPI_Request)*sendCount); 
    int requestIndex = 0;
    // Tag for communication
    int connectivityTag = 7132;
    int coordTag = 7133;
    int eptrTag = 7134;
    int pidTag = 7135;
    int eTypeTag = 7136;

    for (int i = 0; i < world_size; ++i) {
      // Check if there are elements to send
      if (!elementRedistributePattern[i]) {
        continue;
      }
      // Do not send own data
      if (i == world_rank) {
        continue;
      }
      int startLocation = nodeSendCum[i];
      int size = nodeCount[i];
      int sizeElem = elementRedistributePattern[i];
      int startElemLocation = elementSendCum[i];
      // Send connectivity array
      MPI_Isend(&(connectivityPacked[startLocation]), size, MPI_INT, i, \
          connectivityTag, MPI_COMM_WORLD, &(requestListSend[requestIndex]));
      MPI_Isend(&(coordinatesPacked[startLocation*ndim]), ndim*size, MPI_DOUBLE, i, \
          coordTag, MPI_COMM_WORLD, &(requestListSend[requestIndex+1]));
      MPI_Isend(&(eptrPacked[startElemLocation]), sizeElem, MPI_INT, i, \
          eptrTag, MPI_COMM_WORLD, &(requestListSend[requestIndex+2]));
      MPI_Isend(&(pidPacked[startElemLocation]), sizeElem, MPI_INT, i, \
          pidTag, MPI_COMM_WORLD, &(requestListSend[requestIndex+3]));
      MPI_Isend(&(ElementTypePacked[startElemLocation*MAX_ELEMENT_TYPE_SIZE]), \
          sizeElem*MAX_ELEMENT_TYPE_SIZE, MPI_CHAR, i, eTypeTag, MPI_COMM_WORLD, \
          &(requestListSend[requestIndex+4]));
      requestIndex = requestIndex + 5;
    }
    // Verify that number of requests send matches the required
    assert(requestIndex == sendCount);

    // Setup for MPI Receive
    // Create cumulative array to act as recvbuffer pointer
    int *elementRecvCum = (int*)malloc((world_size+1)*sizeof(int));
    int *nodeRecvCum = (int*)malloc((world_size+1)*sizeof(int));
    elementRecvCum[0] = 0;
    nodeRecvCum[0] = 0;
    for (int i = 0; i < world_size; ++i) {
      if (elementRedistributePattern[i] && i != world_rank) {
        recvCount += 1;
      }
      elementRecvCum[i+1] = elementRecvCum[i]+elementReceivePattern[i];
      nodeRecvCum[i+1] = nodeRecvCum[i]+nodeReceivePattern[i];
    }
    // Multiply recvCount by number of messages to recv per process
    recvCount *= 5;
    int nnodesRecv = nodeRecvCum[world_size];
    int nelementsRecv = elementRecvCum[world_size];

    // Create arrays to receive
    int *connectivityRecv = (int*)malloc(nnodesRecv*sizeof(int));
    double *coordinatesRecv = (double*)malloc(nnodesRecv*sizeof(double)*ndim);
    int *eptrRecv = (int*)malloc((nelementsRecv+1)*sizeof(int));
    eptrRecv[0] = 0;
    int *pidRecv = (int*)malloc(nelementsRecv*sizeof(int));
    char *ElementTypeRecv = (char*)malloc(nelementsRecv*\
        MAX_ELEMENT_TYPE_SIZE*sizeof(char));

    MPI_Request* requestListRecv = (MPI_Request*)malloc(sizeof(MPI_Request)*recvCount); 
    requestIndex = 0;
    // Receive arrays
    for (int i = 0; i < world_size; ++i) {
      // If no elements to receive from process i 
      if (!elementReceivePattern[i]) {
        continue;
      }
      int startLocation = nodeRecvCum[i];
      int size = nodeReceivePattern[i];
      int sizeElem = elementReceivePattern[i];
      int startElemLocation = elementRecvCum[i];
      // copy own data
      if (i == world_rank) {
        int startLocSend = nodeSendCum[i];
        int startLocSendElem = elementSendCum[i];
        memcpy(&(connectivityRecv[startLocation]), \
            &(connectivityPacked[startLocSend]), sizeof(int)*size);
        memcpy(&(coordinatesRecv[startLocation*ndim]), \
            &(coordinatesPacked[startLocSend*ndim]), ndim*sizeof(double)*size);
        memcpy(&(eptrRecv[startElemLocation+1]), \
            &(eptrPacked[startLocSendElem]), sizeof(int)*sizeElem);
        memcpy(&(pidRecv[startElemLocation]), \
            &(pidPacked[startLocSendElem]), sizeof(int)*sizeElem);
        memcpy(&(ElementTypeRecv[startElemLocation*MAX_ELEMENT_TYPE_SIZE]), \
            &(ElementTypePacked[startLocSendElem*MAX_ELEMENT_TYPE_SIZE]), \
            sizeof(char)*sizeElem*MAX_ELEMENT_TYPE_SIZE);
        continue;
      }
      // Recv connectivity array from i
      MPI_Irecv(&(connectivityRecv[startLocation]), size, MPI_INT, i, \
          connectivityTag, MPI_COMM_WORLD, &(requestListRecv[requestIndex]));
      MPI_Irecv(&(coordinatesRecv[startLocation*ndim]), ndim*size, MPI_DOUBLE, i, \
          coordTag, MPI_COMM_WORLD, &(requestListRecv[requestIndex+1]));
      MPI_Irecv(&(eptrRecv[startElemLocation+1]), sizeElem, MPI_INT, i, \
          eptrTag, MPI_COMM_WORLD, &(requestListRecv[requestIndex+2]));
      MPI_Irecv(&(pidRecv[startElemLocation]), sizeElem, MPI_INT, i, \
          pidTag, MPI_COMM_WORLD, &(requestListRecv[requestIndex+3]));
      MPI_Irecv(&(ElementTypeRecv[startElemLocation*MAX_ELEMENT_TYPE_SIZE]), \
          sizeElem*MAX_ELEMENT_TYPE_SIZE, MPI_CHAR, i, eTypeTag, \
          MPI_COMM_WORLD, &(requestListRecv[requestIndex+4]));
      requestIndex = requestIndex + 5;
    }
    // Verify that number of requests send matches the required
    assert(requestIndex == recvCount);

    // Wait for all send to be over
    MPI_Status status;
    for (int i = 0; i < sendCount; ++i) {
      MPI_Wait(&(requestListSend[i]), &status);
    }

    if (debug && 1==0) {
      printf("Connectivity Array : \n");
      if (world_rank == 0) {
        printf("Send connectivity array from 0 to 1\n");
        for(int i = 0; i < 10; ++i) {
          printf("%d\t", connectivityPacked[nodeSendCum[1]+i]);
        }
        printf("\n");
        for(int i = 1; i < 11; ++i) {
          printf("%d\t", connectivityPacked[nodeSendCum[2]-i]);
        }
        printf("\n");
      }
    }
    if (debug && 1 == 0) {
      printf("Coord Array : \n");
      if (world_rank == 0) {
        printf("Send coord array from 0 to 1\n");
        for(int i = 0; i < 10; ++i) {
          printf("%.4f\t", coordinatesPacked[ndim*nodeSendCum[1]+i]);
        }
        printf("\n");
        for(int i = 1; i < 11; ++i) {
          printf("%.4f\t", coordinatesPacked[ndim*nodeSendCum[2]-i]);
        }
        printf("\n");
      }
    }
    if (debug && 1 == 0) {
      printf("Eptr Array : \n");
      if (world_rank == 0) {
        printf("Send eptr array from 0 to 1\n");
        for(int i = 0; i < 10; ++i) {
          printf("%d\t", eptrPacked[elementSendCum[1]+i]);
        }
        printf("\n");
        for(int i = 1; i < 11; ++i) {
          printf("%d\t", eptrPacked[elementSendCum[2]-i]);
        }
        printf("\n");
      }
    }
    if (debug && 1 == 0) {
      printf("PID Array : \n");
      if (world_rank == 0) {
        printf("Send pid array from 0 to 1\n");
        for(int i = 0; i < 10; ++i) {
          printf("%d\t", pidPacked[elementSendCum[1]+i]);
        }
        printf("\n");
        for(int i = 1; i < 11; ++i) {
          printf("%d\t", pidPacked[elementSendCum[2]-i]);
        }
        printf("\n");
      }
    }
    if (debug && 1 == 0) {
      printf("Element Type Array : \n");
      if (world_rank == 0) {
        printf("Send element type array from 0 to 1\n");
        for(int i = 0; i < 10; ++i) {
          printf("%s\t", &(ElementTypePacked[(elementSendCum[1]+i)*MAX_ELEMENT_TYPE_SIZE]));
        }
        printf("\n");
        for(int i = 1; i < 11; ++i) {
          printf("%s\t", &(ElementTypePacked[(elementSendCum[2]-i)*MAX_ELEMENT_TYPE_SIZE]));
        }
        printf("\n");
      }
    }

    free(connectivityPacked);
    free(coordinatesPacked);
    free(eptrPacked);
    free(pidPacked);
    free(ElementTypePacked);
    free(requestListSend);

    free(elementRedistributePattern);
    free(nodeRedistributePattern);
    free(elementSendCum);
    free(nodeSendCum);
    free(elementCount);
    free(nodeCount);

    // Wait for all recv to be over
    for (int i = 0; i < recvCount; ++i) {
      MPI_Wait(&(requestListRecv[i]), &status);
    }
    if (debug && 1==0) {
      printf("Connectivity Array : \n");
      if (world_rank == 1) {
        printf("Recv connectivity array from 1 to 0\n");
        for(int i = 0; i < 10; ++i) {
          printf("%d\t", connectivityRecv[nodeRecvCum[0]+i]);
        }
        printf("\n");
        for(int i = 1; i < 11; ++i) {
          printf("%d\t", connectivityRecv[nodeRecvCum[1]-i]);
        }
        printf("\n");
      }
    }
    if (debug && 1 == 0) {
      printf("Coord Array : \n");
      if (world_rank == 1) {
        printf("Recv coord array from 1 to 0\n");
        for(int i = 0; i < 10; ++i) {
          printf("%.4f\t", coordinatesRecv[ndim*nodeRecvCum[0]+i]);
        }
        printf("\n");
        for(int i = 1; i < 11; ++i) {
          printf("%.4f\t", coordinatesRecv[ndim*nodeRecvCum[1]-i]);
        }
        printf("\n");
      }
    }
    if (debug && 1==0) {
      printf("Eptr Array : \n");
      if (world_rank == 1) {
        printf("Recv eptr array from 1 to 0\n");
        for(int i = 0; i < 10; ++i) {
          printf("%d\t", eptrRecv[elementRecvCum[0]+i+1]);
        }
        printf("\n");
        for(int i = 1; i < 11; ++i) {
          printf("%d\t", eptrRecv[elementRecvCum[1]-i+1]);
        }
        printf("\n");
      }
    }
    if (debug && 1==0) {
      printf("PID Array : \n");
      if (world_rank == 1) {
        printf("Recv pid array from 1 to 0\n");
        for(int i = 0; i < 10; ++i) {
          printf("%d\t", pidRecv[elementRecvCum[0]+i]);
        }
        printf("\n");
        for(int i = 1; i < 11; ++i) {
          printf("%d\t", pidRecv[elementRecvCum[1]-i]);
        }
        printf("\n");
      }
    }
    if (debug && 1 == 0) {
      printf("Element Type Array : \n");
      if (world_rank == 0) {
        printf("Recv element type array from 0 to 1\n");
        for(int i = 0; i < 10; ++i) {
          printf("%s\t", &(ElementTypeRecv[(elementRecvCum[1]+i)*MAX_ELEMENT_TYPE_SIZE]));
        }
        printf("\n");
        for(int i = 1; i < 11; ++i) {
          printf("%s\t", &(ElementTypeRecv[(elementRecvCum[2]-i)*MAX_ELEMENT_TYPE_SIZE]));
        }
        printf("\n");
      }
    }
    free(elementReceivePattern);
    free(nodeReceivePattern);
    free(requestListRecv);

    // Free Recv arrays
    free(elementRecvCum);
    free(nodeRecvCum);

    // Free and swap original arrays
    free(connectivity);
    connectivity = connectivityRecv;
    free(coordinates);
    coordinates = coordinatesRecv;
    free(pid);
    pid = pidRecv;
    
    // Find cumulative sum on eptr
    for (int i = 0; i < nelementsRecv; ++i) {
      eptrRecv[i+1] += eptrRecv[i];
    }
    free(eptr);
    eptr = eptrRecv;

    // Unpack and create Element type
    for (int i = 0; i < nelements; ++i) {
      free(ElementType[i]);
    }
    free(ElementType);
    ElementType = (char**)malloc(nelementsRecv*sizeof(char*));
    for (int i = 0; i < nelementsRecv; ++i) {
      ElementType[i] = (char*)malloc(MAX_ELEMENT_TYPE_SIZE*sizeof(char));
      memcpy(ElementType[i], &(ElementTypeRecv[i*MAX_ELEMENT_TYPE_SIZE]), MAX_ELEMENT_TYPE_SIZE*sizeof(char));
    }
    free(ElementTypeRecv);

    nnodes = eptr[nelementsRecv];
    nelements = nelementsRecv;

    // Reorder local connectivity
    updateConnectivityGlobalToLocal();
    // assert(nelements == 0);
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
