#include "FemTech.h"

#include <assert.h>
#include <stdlib.h>

// Variables to keep store the communication patterns between processes
int sendProcessCount;
int *sendProcessID;
int *sendNeighbourCount;
int *sendNeighbourCountCum;
int *sendNodeIndex;
double *sendNodeDisplacement;
double *recvNodeDisplacement;

int *globalNodeID;

void updateConnectivityGlobalToLocal(void);
void createNodalCommunicationPattern(void);

void PartitionMesh() {

  // Options for ParMETIS
  idx_t options[3];
  options[0] = 1;  // timing info printed
  options[1] = 0;  // timing info suppressed
  options[2] = 15; // seed for random number generator

  // Number of partitions (one for each process)
  idx_t npartis = world_size; // number of parts equals number of processes

  // Strange weight arrays needed by ParMETIS
  idx_t ncon = 1; // number of balance constraints

  // Prepare remaining arguments for ParMETIS
  // idx_t *elmwgt = NULL;
  idx_t *elmwgt = (idx_t*)malloc(nelements*sizeof(idx_t));
  // Assign weights based on element type and material type
  for (int i = 0; i < nelements; ++i) {
    const unsigned int matID = materialID[pid[i]];
    const int matW = materialWeight[matID];
    const int elemW = elemWeight[elemID.at(ElementType[i])];
    elmwgt[i] = matW*elemW;
  }
  idx_t wgtflag = 2;      // we don't use weights
  idx_t numflag = 0;      // we are using C-style arrays
  idx_t ncommonnodes = 2; // number of nodes elements must have in common
  if (nallelements == 2) {
    FILE_LOG(WARNING, "PartitionMesh.cpp: nallelements = 2 so ncommonnodes set to 8 (for "
           "parmetis)");
    ncommonnodes = 8;
  }
  const int tpwgts_size = ncon * npartis;
  real_t *tpwgts = (real_t *)malloc(tpwgts_size * sizeof(real_t));
  for (int i = 0; i < tpwgts_size; i++) {
    tpwgts[i] = (real_t)1.0 / (real_t)npartis;
  }

  real_t *ubvec = (real_t *)malloc(ncon * sizeof(real_t));
  for (int i = 0; i < ncon; i++) {
    ubvec[i] = (real_t)1.05; // weight tolerance (same weight tolerance for
                             // every weight there is)
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
  // Size of this array equals to number/count of mesh elements local to
  // processor This array indicates which processor i-th element belongs to
  idx_t edgecut; // output of ParMETIS_V3_PartMeshKway function
  idx_t *part = (idx_t *)calloc(nelements, sizeof(idx_t));

  MPI_Comm comm;
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);
  const int Result = ParMETIS_V3_PartMeshKway(
      elmdist, eptr, connectivity, elmwgt, &wgtflag, &numflag, &ncon,
      &ncommonnodes, &npartis, tpwgts, ubvec, options, &edgecut, part, &comm);

  // Output of ParMETIS_V3_PartMeshKway call will be used here
  if (Result == METIS_OK) {
    // Variable to count number of elements to move to each MPI process
    // from the current local set of elements
    // variable part the output of parmetis is used to populate this array
    int *elementRedistributePattern = (int *)calloc(world_size, sizeof(int));
    int *nodeRedistributePattern = (int *)calloc(world_size, sizeof(int));
    int processorToMoveElement, cordInElement;
    for (int i = 0; i < nelements; ++i) {
      processorToMoveElement = (int)part[i];
      cordInElement = eptr[i + 1] - eptr[i];
      elementRedistributePattern[processorToMoveElement] += 1;
      nodeRedistributePattern[processorToMoveElement] += cordInElement;
    }
    // Create cumulative array to act as sendbuffer pointer
    int *elementSendCum = (int *)malloc((world_size + 1) * sizeof(int));
    int *nodeSendCum = (int *)malloc((world_size + 1) * sizeof(int));
    elementSendCum[0] = 0;
    nodeSendCum[0] = 0;
    // Variables to keep track of requests
    int sendCount = 0, recvCount = 0;
    for (int i = 0; i < world_size; ++i) {
      if (elementRedistributePattern[i] && i != world_rank) {
        sendCount += 1;
      }
      elementSendCum[i + 1] = elementSendCum[i] + elementRedistributePattern[i];
      nodeSendCum[i + 1] = nodeSendCum[i] + nodeRedistributePattern[i];
    }
    // Multiply by 6 as we are sending 6 seperate lists for each process
    sendCount *= 6;
    int nnodes_process = nodeSendCum[world_size];
    // Verify that all elements are included for distribution
    assert(nelements == elementSendCum[world_size]);

    // Get count to receive from each processor
    int *elementReceivePattern = (int *)malloc(world_size * sizeof(int));
    int *nodeReceivePattern = (int *)malloc(world_size * sizeof(int));
    memcpy(elementReceivePattern, elementRedistributePattern,
           world_size * sizeof(int));
    memcpy(nodeReceivePattern, nodeRedistributePattern,
           world_size * sizeof(int));
    MPI_Alltoall(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, elementReceivePattern, 1,
                 MPI_INT, MPI_COMM_WORLD);
    MPI_Alltoall(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, nodeReceivePattern, 1,
                 MPI_INT, MPI_COMM_WORLD);

    FILE_LOGArrayInt(DEBUGLOGIGNORE, elementRedistributePattern, world_size, \
        "Element Redistribution Pattern");
    FILE_LOGArrayInt(DEBUGLOGIGNORE, elementReceivePattern, world_size, \
        "Element Receive Pattern");
    FILE_LOGArrayInt(DEBUGLOGIGNORE, nodeRedistributePattern, world_size, \
        "Node Redistribution Pattern");
    FILE_LOGArrayInt(DEBUGLOGIGNORE, nodeReceivePattern, world_size, \
        "Node Receive Pattern");

    // Pack data to redistribute mesh
    // Data to send and receive include node coordinates, element connectivity,
    // eptr array, part id array, Element type array
    // Create send buffers

    // Create array to keep track of current count during packing
    int *elementCount = (int *)calloc(world_size, sizeof(int));
    int *nodeCount = (int *)calloc(world_size, sizeof(int));

    // Variables to reorder local data to send
    int *connectivityPacked = (int *)malloc(nnodes_process * sizeof(int));
    double *coordinatesPacked =
        (double *)malloc(nnodes_process * sizeof(double) * ndim);
    int *eptrPacked = (int *)malloc(nelements * sizeof(int));
    int *pidPacked = (int *)malloc(nelements * sizeof(int));
    int *eidPacked = (int *)malloc(nelements * sizeof(int));
    char *ElementTypePacked =
        (char *)malloc(nelements * MAX_ELEMENT_TYPE_SIZE * sizeof(char *));

    // Move data to packed arrays
    for (int i = 0; i < nelements; ++i) {
      // Anil(TODO) This can be optimized to move multiple elements if
      // concecutive elements to move are assigned same destinaation process
      int start = eptr[i];
      int end = eptr[i + 1];
      processorToMoveElement = (int)part[i];
      int count = end - start;
      int locationToCopy = nodeSendCum[processorToMoveElement] +
                           nodeCount[processorToMoveElement];
      int elementLocToCopy = elementSendCum[processorToMoveElement] +
                             elementCount[processorToMoveElement];
      // Copy connectivity
      memcpy(&(connectivityPacked[locationToCopy]), &(connectivity[start]),
             sizeof(int) * count);
      // Copy coordinates
      memcpy(&(coordinatesPacked[locationToCopy * ndim]),
             &(coordinates[start * ndim]), sizeof(double) * count * ndim);
      // Store number of nodes to recreate eptr on destination
      eptrPacked[elementLocToCopy] = count;
      // Copy pid to send
      pidPacked[elementLocToCopy] = pid[i];
      eidPacked[elementLocToCopy] = global_eid[i];
      // Copy element type
      memcpy(&(ElementTypePacked[elementLocToCopy * MAX_ELEMENT_TYPE_SIZE]),
             ElementType[i], sizeof(char) * MAX_ELEMENT_TYPE_SIZE);

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
    MPI_Request *requestListSend =
        (MPI_Request *)malloc(sizeof(MPI_Request) * sendCount);
    int requestIndex = 0;
    // Tag for communication
    int connectivityTag = 7132;
    int coordTag = 7133;
    int eptrTag = 7134;
    int pidTag = 7135;
    int eTypeTag = 7136;
    int eidTag = 7137;

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
      MPI_Isend(&(connectivityPacked[startLocation]), size, MPI_INT, i,
                connectivityTag, MPI_COMM_WORLD,
                &(requestListSend[requestIndex]));
      MPI_Isend(&(coordinatesPacked[startLocation * ndim]), ndim * size,
                MPI_DOUBLE, i, coordTag, MPI_COMM_WORLD,
                &(requestListSend[requestIndex + 1]));
      MPI_Isend(&(eptrPacked[startElemLocation]), sizeElem, MPI_INT, i, eptrTag,
                MPI_COMM_WORLD, &(requestListSend[requestIndex + 2]));
      MPI_Isend(&(pidPacked[startElemLocation]), sizeElem, MPI_INT, i, pidTag,
                MPI_COMM_WORLD, &(requestListSend[requestIndex + 3]));
      MPI_Isend(&(ElementTypePacked[startElemLocation * MAX_ELEMENT_TYPE_SIZE]),
                sizeElem * MAX_ELEMENT_TYPE_SIZE, MPI_CHAR, i, eTypeTag,
                MPI_COMM_WORLD, &(requestListSend[requestIndex + 4]));
      MPI_Isend(&(eidPacked[startElemLocation]), sizeElem, MPI_INT, i, eidTag,
                MPI_COMM_WORLD, &(requestListSend[requestIndex + 5]));
      requestIndex = requestIndex + 6;
    }
    // Verify that number of requests send matches the required
    assert(requestIndex == sendCount);

    // Setup for MPI Receive
    // Create cumulative array to act as recvbuffer pointer
    int *elementRecvCum = (int *)malloc((world_size + 1) * sizeof(int));
    int *nodeRecvCum = (int *)malloc((world_size + 1) * sizeof(int));
    elementRecvCum[0] = 0;
    nodeRecvCum[0] = 0;
    for (int i = 0; i < world_size; ++i) {
      if (elementReceivePattern[i] && i != world_rank) {
        recvCount += 1;
      }
      elementRecvCum[i + 1] = elementRecvCum[i] + elementReceivePattern[i];
      nodeRecvCum[i + 1] = nodeRecvCum[i] + nodeReceivePattern[i];
    }
    // Multiply recvCount by number of messages to recv per process
    recvCount *= 6;
    int nnodesRecv = nodeRecvCum[world_size];
    int nelementsRecv = elementRecvCum[world_size];

    // Create arrays to receive
    int *connectivityRecv = (int *)malloc(nnodesRecv * sizeof(int));
    double *coordinatesRecv =
        (double *)malloc(nnodesRecv * sizeof(double) * ndim);
    int *eptrRecv = (int *)malloc((nelementsRecv + 1) * sizeof(int));
    eptrRecv[0] = 0;
    int *pidRecv = (int *)malloc(nelementsRecv * sizeof(int));
    int *eidRecv = (int *)malloc(nelementsRecv * sizeof(int));
    char *ElementTypeRecv =
        (char *)malloc(nelementsRecv * MAX_ELEMENT_TYPE_SIZE * sizeof(char));

    MPI_Request *requestListRecv =
        (MPI_Request *)malloc(sizeof(MPI_Request) * recvCount);
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
        memcpy(&(connectivityRecv[startLocation]),
               &(connectivityPacked[startLocSend]), sizeof(int) * size);
        memcpy(&(coordinatesRecv[startLocation * ndim]),
               &(coordinatesPacked[startLocSend * ndim]),
               ndim * sizeof(double) * size);
        memcpy(&(eptrRecv[startElemLocation + 1]),
               &(eptrPacked[startLocSendElem]), sizeof(int) * sizeElem);
        memcpy(&(pidRecv[startElemLocation]), &(pidPacked[startLocSendElem]),
               sizeof(int) * sizeElem);
        memcpy(&(ElementTypeRecv[startElemLocation * MAX_ELEMENT_TYPE_SIZE]),
               &(ElementTypePacked[startLocSendElem * MAX_ELEMENT_TYPE_SIZE]),
               sizeof(char) * sizeElem * MAX_ELEMENT_TYPE_SIZE);
        memcpy(&(eidRecv[startElemLocation]), &(eidPacked[startLocSendElem]),
               sizeof(int) * sizeElem);
        continue;
      }
      // Recv connectivity array from i
      MPI_Irecv(&(connectivityRecv[startLocation]), size, MPI_INT, i,
                connectivityTag, MPI_COMM_WORLD,
                &(requestListRecv[requestIndex]));
      MPI_Irecv(&(coordinatesRecv[startLocation * ndim]), ndim * size,
                MPI_DOUBLE, i, coordTag, MPI_COMM_WORLD,
                &(requestListRecv[requestIndex + 1]));
      MPI_Irecv(&(eptrRecv[startElemLocation + 1]), sizeElem, MPI_INT, i,
                eptrTag, MPI_COMM_WORLD, &(requestListRecv[requestIndex + 2]));
      MPI_Irecv(&(pidRecv[startElemLocation]), sizeElem, MPI_INT, i, pidTag,
                MPI_COMM_WORLD, &(requestListRecv[requestIndex + 3]));
      MPI_Irecv(&(ElementTypeRecv[startElemLocation * MAX_ELEMENT_TYPE_SIZE]),
                sizeElem * MAX_ELEMENT_TYPE_SIZE, MPI_CHAR, i, eTypeTag,
                MPI_COMM_WORLD, &(requestListRecv[requestIndex + 4]));
      MPI_Irecv(&(eidRecv[startElemLocation]), sizeElem, MPI_INT, i, eidTag,
                MPI_COMM_WORLD, &(requestListRecv[requestIndex + 5]));
      requestIndex = requestIndex + 6;
    }
    // Verify that number of requests send matches the required
    assert(requestIndex == recvCount);

    // Wait for all send to be over
    MPI_Status status;
    for (int i = 0; i < sendCount; ++i) {
      MPI_Wait(&(requestListSend[i]), &status);
    }

    FILE_LOGArrayIntMaster(DEBUGLOGIGNORE, &connectivityPacked[nodeSendCum[1]], \
        10, "Send connectivity array from 0 to 1");
    FILE_LOGArrayIntMaster(DEBUGLOGIGNORE, &connectivityPacked[nodeSendCum[2]-10], \
        10, "Send connectivity array from 0 to 1");
    FILE_LOGArrayMaster(DEBUGLOGIGNORE, &coordinatesPacked[ndim*nodeSendCum[1]], \
        10, "Send coord array from 0 to 1");
    FILE_LOGArrayMaster(DEBUGLOGIGNORE, &coordinatesPacked[ndim*nodeSendCum[2]-10], \
        10, "Send coord array from 0 to 1");
    FILE_LOGArrayIntMaster(DEBUGLOGIGNORE, &eptrPacked[elementSendCum[1]], \
        10, "Send eptr array from 0 to 1");
    FILE_LOGArrayIntMaster(DEBUGLOGIGNORE, &eptrPacked[elementSendCum[2]-10], \
        10, "Send eptr array from 0 to 1");
    FILE_LOGArrayIntMaster(DEBUGLOGIGNORE, &pidPacked[elementSendCum[1]], \
        10, "Send pid array from 0 to 1");
    FILE_LOGArrayIntMaster(DEBUGLOGIGNORE, &pidPacked[elementSendCum[2]-10], \
        10, "Send pid array from 0 to 1");

    free(connectivityPacked);
    free(coordinatesPacked);
    free(eptrPacked);
    free(pidPacked);
    free(eidPacked);
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

    if (world_rank == 1) {
      FILE_LOGArrayIntSingle(DEBUGLOGIGNORE, &connectivityRecv[nodeRecvCum[0]], \
          10, "Recv connectivity array from 1 to 0");
      FILE_LOGArrayIntSingle(DEBUGLOGIGNORE, &connectivityRecv[nodeRecvCum[1]-10], \
          10, "Recv connectivity array from 1 to 0");
      FILE_LOGArraySingle(DEBUGLOGIGNORE, &coordinatesRecv[ndim*nodeRecvCum[0]], \
          10, "Recv coord array from 1 to 0");
      FILE_LOGArraySingle(DEBUGLOGIGNORE, &coordinatesRecv[ndim*nodeRecvCum[1]-10], \
          10, "Recv coord array from 1 to 0");
      FILE_LOGArrayIntSingle(DEBUGLOGIGNORE, &eptrRecv[elementRecvCum[0]+1], \
          10, "Recv eptr array from 1 to 0");
      FILE_LOGArrayIntSingle(DEBUGLOGIGNORE, &eptrRecv[elementRecvCum[1]-9], \
          10, "Recv eptr array from 1 to 0");
      FILE_LOGArrayIntSingle(DEBUGLOGIGNORE, &pidRecv[elementRecvCum[0]], \
          10, "Recv pid array from 1 to 0");
      FILE_LOGArrayIntSingle(DEBUGLOGIGNORE, &pidRecv[elementRecvCum[1]-10], \
          10, "Recv pid array from 1 to 0");
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
    free(global_eid);
    global_eid = eidRecv;

    // Find cumulative sum on eptr
    for (int i = 0; i < nelementsRecv; ++i) {
      eptrRecv[i + 1] += eptrRecv[i];
    }
    free(eptr);
    eptr = eptrRecv;

    // Unpack and create Element type
    for (int i = 0; i < nelements; ++i) {
      free(ElementType[i]);
    }
    free(ElementType);
    ElementType = (char **)malloc(nelementsRecv * sizeof(char *));
    for (int i = 0; i < nelementsRecv; ++i) {
      ElementType[i] = (char *)malloc(MAX_ELEMENT_TYPE_SIZE * sizeof(char));
      memcpy(ElementType[i], &(ElementTypeRecv[i * MAX_ELEMENT_TYPE_SIZE]),
             MAX_ELEMENT_TYPE_SIZE * sizeof(char));
    }
    free(ElementTypeRecv);

    nNodes = eptr[nelementsRecv];
    nelements = nelementsRecv;

    // Create nodal communcication pattern
    createNodalCommunicationPattern();

    // Reorder local connectivity
    updateConnectivityGlobalToLocal();
    nDOF = nNodes * ndim;

    FILE_LOG(INFO, "Number of nodes : %d", nNodes);
    FILE_LOG(INFO, "Number of DOFs : %d", nDOF);
  } else {
    FILE_LOG_SINGLE(ERROR, "ParMETIS returned error code %d", Result);
    exit(EXIT_FAILURE);
  }
  free(ubvec);
  free(tpwgts);
  free(elmdist);
  free(part);
  free(elmwgt);
}
//-------------------------------------------------------------------------------------------
int compare(const void *a, const void *b) { return (*(int *)a - *(int *)b); }
//-------------------------------------------------------------------------------------------
int unique(int *arr, int n) {
  int *temp;
  int j = 0;
  temp = (int *)calloc(n, sizeof(int));
  for (int i = 0; i < n - 1; i++) {
    if (arr[i] != arr[i + 1]) {
      temp[j++] = arr[i];
    }
  }
  temp[j++] = arr[n - 1];
  for (int i = 0; i < j; i++) {
    arr[i] = temp[i];
  }
  free(temp);
  return j;
}
int coordinateFromGlobalID(int *array, int nodeID, int size, double* coord) {
  int *position = (int*)bsearch(&nodeID, array, size, sizeof(int), compare);
  if (position) {
    int location = position-array;
    // printf("Value at position : %d or %d\n", *position, array[location]);
    // printf("NodeID : %d (location = %d), found on rank %d\n", nodeID, location, world_rank);
    location *= 3;
    coord[0] = coordinates[location]; 
    coord[1] = coordinates[location+1]; 
    coord[2] = coordinates[location+2]; 
    // printf("Coordinates : %f, %f, %f\n", coord[0], coord[1], coord[2]);
    return 1;
  } else {
    return 0;
  }
}
//-------------------------------------------------------------------------------------------
void updateConnectivityGlobalToLocal(void) {
  int totalSize = eptr[nelements];
  int *newConnectivity = (int *)malloc(totalSize * sizeof(int));
  int *sorted = (int *)malloc(totalSize * sizeof(int));
  memcpy(sorted, connectivity, totalSize * sizeof(int));
  qsort(sorted, totalSize, sizeof(int), compare);
  nNodes = unique(sorted, totalSize);
  for (int i = 0; i < totalSize; ++i) {
    int j;
    for (j = 0; j < nNodes; ++j) {
      if (sorted[j] == connectivity[i]) {
        break;
      }
    }
    newConnectivity[i] = j;
  }

  // Reoder co-ordinates
  double *newCoordinates = (double *)malloc(ndim * nNodes * sizeof(double));
  for (int j = 0; j < nNodes; ++j) {
    int i;
    for (i = 0; i < totalSize; ++i) {
      if (sorted[j] == connectivity[i]) {
        break;
      }
    }
    memcpy(newCoordinates + ndim * j, coordinates + ndim * i,
           ndim * sizeof(double));
  }
  // Global to local map of sendNodeIndex Array
  totalSize = sendNeighbourCountCum[sendProcessCount];
  for (int i = 0; i < totalSize; ++i) {
    int j;
    for (j = 0; j < nNodes; ++j) {
      if (sorted[j] == sendNodeIndex[i]) {
        break;
      }
    }
    sendNodeIndex[i] = j;
  }
  // Copy sorted array to globalNodeID for pre and post processing
  globalNodeID = (int*)malloc(nNodes*sizeof(int));
  for (int i = 0; i < nNodes; ++i) {
    globalNodeID[i] = sorted[i];
  }
  free(connectivity);
  connectivity = newConnectivity;
  free(coordinates);
  coordinates = newCoordinates;
  free(sorted);
}
// Function to compute intersection of two sorted arrays
int getIntersection(int *list1, int *list2, int size1, int size2,
                    int **intersection) {
  // Create a temperory array of size of smallest of the two arrays
  int maxSize = (size1 < size2) ? size1 : size2;
  int *temp = (int *)malloc(maxSize * sizeof(int));
  int i = 0, j = 0, k = 0;

  while (i < size1 && j < size2) {
    if (list1[i] < list2[j]) {
      i = i + 1;
    } else {
      if (list2[j] < list1[i]) {
        j = j + 1;
      } else { // Both list are equal
        temp[k] = list2[j];
        i = i + 1;
        j = j + 1;
        k = k + 1;
      }
    }
  }

  int intSize = k; // Size of the intersection
  (*intersection) = (int *)malloc(intSize * sizeof(int));
  memcpy(*intersection, temp, intSize * sizeof(int));
  free(temp);
  return intSize;
}
void createNodalCommunicationPattern(void) {
  // Creates the nodal communication pattern between processes
  // 1 . Identify the neighbours or ghost node of each element after the mesh is
  // partitoned.
  // 2. Find the boundary and interior nodes.
  // 3. For boundary nodes identify the process ID of neighbour ghost elements
  // 4. Send element ID to neighbour process
  // 5. Create list of node list of received element ID
  // 6. Send back the node list and node ptr of ghost elements
  // 7. Create a list of common node and the processID to which it belong
  // 8. Finally populate the data structure that can be used for communication
  // every time step.

  // Create dual of mesh to figure out ghost elements of the local mesh
  MPI_Comm comm;
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);
  idx_t *xadj;
  idx_t *adjncy;
  idx_t ncommonnodes = 1; // number of nodes elements must have in common
  idx_t numflag = 0;      // we are using C-style arrays
  // Create the elelemnt distribution array
  int *elementCount = (int *)malloc(world_size * sizeof(int));
  MPI_Allgather(&nelements, 1, MPI_INT, elementCount, 1, MPI_INT, comm);

  idx_t *elmdist = (idx_t *)malloc((world_size + 1) * sizeof(idx_t));
  elmdist[0] = 0;
  for (int i = 0; i < world_size; ++i) {
    elmdist[i + 1] = elmdist[i] + elementCount[i];
  }
#ifdef DEBUG
  for (int i = 0; i < world_size; ++i) {
    FILE_LOG(DEBUGLOGIGNORE, "%d, %d", elementCount[i], elmdist[i + 1]);
  }
#endif //DEBUG
  free(elementCount);
  const int Result =
      ParMETIS_V3_Mesh2Dual(elmdist, eptr, connectivity, &numflag,
                            &ncommonnodes, &xadj, &adjncy, &comm);
  // Check the success of Mesh2Dual
  assert(Result == METIS_OK);

  // Partition the neighbours to process they belong to
  const int totalGhosts = xadj[nelements];
  int *processID_ghost = (int *)malloc(totalGhosts * sizeof(int));

  int *ghostCountProcess = (int *)calloc(world_size, sizeof(int));
  int *ghostCountProcessCum = (int *)malloc((world_size + 1) * sizeof(int));
  int totalGhostCount = 0;

  for (int i = 0; i < totalGhosts; ++i) {
    int elementID = adjncy[i];
    int j;
    for (j = 0; j < world_size; ++j) {
      if (elementID < elmdist[j + 1]) {
        processID_ghost[i] = j;
        if (j != world_rank) {
          totalGhostCount += 1;
          ghostCountProcess[j] += 1;
        }
        break;
      }
    }
    // Confirm that a processID has been assigned
    assert(j != world_size);
  }

  ghostCountProcessCum[0] = 0;
  for (int i = 0; i < world_size; ++i) {
    ghostCountProcessCum[i + 1] =
        ghostCountProcessCum[i] + ghostCountProcess[i];
  }

  assert(ghostCountProcessCum[world_size] == totalGhostCount);

  // Classify elements as interior or boundary based on process id of ghost
  // If all the neighbours of a cell are in the same process its classified as
  // an interior element
  // This can be used to skip all interior elements quickly
  int *boundaryElement = (int *)calloc(nelements, sizeof(int));
  for (int i = 0; i < nelements; ++i) {
    int end = xadj[i + 1];
    for (int j = xadj[i]; j < end; ++j) {
      if (processID_ghost[j] != world_rank) {
        boundaryElement[i] = 1;
        break;
      }
    }
  }

  // Create array to store elementID for request from other process
  int *elemID_request = (int *)malloc(totalGhostCount * sizeof(int));
  int *elemCurrentLoc = (int *)calloc(world_size, sizeof(int));

  // Node list count of elements that share with ghost nodes
  // This count will be used in the last step to create local node list to
  // compare with ghost node list from each process
  int *nodeListCountProcess = (int *)calloc(world_size, sizeof(int));
  for (int j = 0; j < nelements; ++j) {
    if (boundaryElement[j]) {
      int localElementNodeSize = eptr[j + 1] - eptr[j];
      int end = xadj[j + 1];
      for (int i = xadj[j]; i < end; ++i) {
        int elementID = adjncy[i];
        int processID = processID_ghost[i];
        if (processID != world_rank) {
          int location =
              ghostCountProcessCum[processID] + elemCurrentLoc[processID];
          elemID_request[location] = elementID;
          elemCurrentLoc[processID] += 1;
          nodeListCountProcess[processID] += localElementNodeSize;
        }
      }
    }
  }
  // Check if all ghosts are added
  for (int i = 0; i < world_size; ++i) {
    assert(elemCurrentLoc[i] == ghostCountProcess[i]);
  }
  free(elemCurrentLoc);

  // Create unique set of elementIDs to request by removing duplicates in
  // elemID_request
  // Array ghostCountProcess is updated to store the count of unique element ids
  int sendCount = 0;
  for (int i = 0; i < world_size; ++i) {
    if (ghostCountProcess[i]) {
      sendCount += 1;
#ifdef DEBUG
      if (debug && 1 == 0) {
        printf("DEBUG (%d) : Before Make Unique for %d\n", world_rank, i);
        for (int j = ghostCountProcessCum[i]; j < ghostCountProcessCum[i + 1];
             ++j) {
          printf("%d\t", elemID_request[j]);
        }
        printf("\n");
      }
#endif //DEBUG
      qsort(&(elemID_request[ghostCountProcessCum[i]]), ghostCountProcess[i],
            sizeof(int), compare);
      ghostCountProcess[i] = unique(&(elemID_request[ghostCountProcessCum[i]]),
                                    ghostCountProcess[i]);
#ifdef DEBUG
      if (debug && 1 == 0) {
        printf("DEBUG (%d) : After Make Unique for %d\n", world_rank, i);
        for (int j = 0; j < ghostCountProcess[i]; ++j) {
          printf("%d\t", elemID_request[ghostCountProcessCum[i] + j]);
        }
        printf("\n");
      }
#endif //DEBUG
    }
  }
  // Prepare arrays to recv request for ghost
  // Get count to receive from each processor
  int *elemID_recvCount = (int *)malloc(world_size * sizeof(int));
  memcpy(elemID_recvCount, ghostCountProcess, world_size * sizeof(int));
  MPI_Alltoall(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, elemID_recvCount, 1, MPI_INT,
               MPI_COMM_WORLD);

  FILE_LOGArrayInt(DEBUGLOGIGNORE, ghostCountProcess, world_size, \
      "Element ID Request Pattern");
  FILE_LOGArrayInt(DEBUGLOGIGNORE, elemID_recvCount, world_size, \
      "Element ID Receive Pattern");

  // Create the required send and receive patterns
  int *elemID_recvCountCum = (int *)malloc((world_size + 1) * sizeof(int));
  elemID_recvCountCum[0] = 0;

  int recvCount = 0;
  for (int i = 0; i < world_size; ++i) {
    elemID_recvCountCum[i + 1] = elemID_recvCountCum[i] + elemID_recvCount[i];
    if (elemID_recvCount[i]) {
      recvCount += 1;
    }
  }
  int totalElemID_ToRecv = elemID_recvCountCum[world_size];
  int *elemID_recv = (int *)malloc(totalElemID_ToRecv * sizeof(int));

  // factor of 2 used in below malloc required to transfer nodeID and nodePtr
  MPI_Request *requestListSend =
      (MPI_Request *)malloc(sizeof(MPI_Request) * 2 * sendCount);
  MPI_Request *requestListRecv =
      (MPI_Request *)malloc(sizeof(MPI_Request) * 2 * recvCount);
  int requestIndex = 0;
  int elemIDTag = 923;

  // Post all send requests
  for (int i = 0; i < world_size; ++i) {
    if (ghostCountProcess[i]) {
      MPI_Isend(&(elemID_request[ghostCountProcessCum[i]]),
                ghostCountProcess[i], MPI_INT, i, elemIDTag, MPI_COMM_WORLD,
                &(requestListSend[requestIndex]));
      requestIndex += 1;
    }
  }
  // Check if sends tally
  assert(requestIndex == sendCount);

  // Post all recv requests
  requestIndex = 0;
  for (int i = 0; i < world_size; ++i) {
    if (elemID_recvCount[i]) {
      // Recv elemt id list array from i
      MPI_Irecv(&(elemID_recv[elemID_recvCountCum[i]]), elemID_recvCount[i],
                MPI_INT, i, elemIDTag, MPI_COMM_WORLD,
                &(requestListRecv[requestIndex]));
      requestIndex += 1;
    }
  }
  // Check if recvs tally
  assert(requestIndex == recvCount);

  // Wait for all sends and recvs to be over
  MPI_Status status;
  for (int i = 0; i < sendCount; ++i) {
    MPI_Wait(&(requestListSend[i]), &status);
  }
  for (int i = 0; i < recvCount; ++i) {
    MPI_Wait(&(requestListRecv[i]), &status);
  }

#ifdef DEBUG
  if (debug && 1 == 0) {
    printf("DEBUG(%d) : Element id requests received \n", world_rank);
    for (int i = 0; i < world_size; ++i) {
      printf("From %d\n", i);
      for (int j = elemID_recvCountCum[i]; j < elemID_recvCountCum[i + 1];
           ++j) {
        printf("%d\t", elemID_recv[j]);
      }
      printf("\n");
    }
  }
#endif //DEBUG

  // Prepare to send nodelist of requested elements
  int totalElementRequest = elemID_recvCountCum[world_size];
  int *nodePtr = (int *)malloc((totalElementRequest + 1) * sizeof(int));
  int *nodeCount = (int *)malloc(totalElementRequest * sizeof(int));
  nodePtr[0] = 0;
  for (int i = 0; i < totalElementRequest; ++i) {
    // Convert global element ID to local element ID
    int elementID = elemID_recv[i] - elmdist[world_rank];
    nodeCount[i] = eptr[elementID + 1] - eptr[elementID];
    nodePtr[i + 1] = nodePtr[i] + nodeCount[i];
  }

  int totalNodes = nodePtr[totalElementRequest];
  // Create node list
  int *nodeList = (int *)malloc(totalNodes * sizeof(int));
  for (int i = 0; i < totalElementRequest; ++i) {
    // Convert global element ID to local element ID
    int elementID = elemID_recv[i] - elmdist[world_rank];
    int nodesInElement = eptr[elementID + 1] - eptr[elementID];
    memcpy(&(nodeList[nodePtr[i]]), &(connectivity[eptr[elementID]]),
           sizeof(int) * nodesInElement);
  }

  free(elmdist);
  free(elemID_recv);

  FILE_LOGArrayInt(DEBUGLOGIGNORE, nodeList, totalNodes, "Node id populated");

  // Send nodesIDs associated with corresponding elemntID to requested process
  // First exchange number of nodeIDs to send between process
  // Create array for the same
  int *nodeID_countSend = (int *)calloc(world_size, sizeof(int));
  for (int i = 0; i < world_size; ++i) {
    if (elemID_recvCount[i]) {
      int start = elemID_recvCountCum[i];
      int end = elemID_recvCountCum[i + 1];
      nodeID_countSend[i] = nodePtr[end] - nodePtr[start];
    }
  }
  int *nodeID_countRecv = (int *)malloc(world_size * sizeof(int));
  int *nodeID_countRecvCum = (int *)malloc((world_size + 1) * sizeof(int));
  memcpy(nodeID_countRecv, nodeID_countSend, world_size * sizeof(int));
  MPI_Alltoall(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, nodeID_countRecv, 1, MPI_INT,
               MPI_COMM_WORLD);

  FILE_LOGArrayInt(DEBUGLOGIGNORE, nodeID_countSend, world_size, "Node ID Send Pattern");
  FILE_LOGArrayInt(DEBUGLOGIGNORE, nodeID_countRecv, world_size, "Node ID Recv Pattern");

  // Create storage space to receive node ID
  nodeID_countRecvCum[0] = 0;
  // Update ghostCountProcessCum
  int *uniqueElementCountCum = (int *)malloc((world_size + 1) * sizeof(int));
  uniqueElementCountCum[0] = 0;
  for (int i = 0; i < world_size; ++i) {
    nodeID_countRecvCum[i + 1] = nodeID_countRecvCum[i] + nodeID_countRecv[i];
    uniqueElementCountCum[i + 1] =
        uniqueElementCountCum[i] + ghostCountProcess[i];
  }
  int nodePtrSize = uniqueElementCountCum[world_size];
  int totalNodeIDToRecv = nodeID_countRecvCum[world_size];

  int *nodeListRecv = (int *)malloc(totalNodeIDToRecv * sizeof(int));
  int *nodePtrRecv = (int *)malloc((nodePtrSize + 1) * sizeof(int));
  nodePtrRecv[0] = 0;

  // Send Node ID to all requested process
  requestIndex = 0;
  int nodeIDTag = 924;
  int nodeID_PtrTag = 925;

  for (int i = 0; i < world_size; ++i) {
    if (elemID_recvCount[i]) {
      int start = elemID_recvCountCum[i];
      int location = nodePtr[start];
      int nodeIDCount = nodeID_countSend[i];
      MPI_Isend(&(nodeList[location]), nodeIDCount, MPI_INT, i, nodeIDTag,
                MPI_COMM_WORLD, &(requestListRecv[requestIndex]));
      MPI_Isend(&(nodeCount[start]), elemID_recvCount[i], MPI_INT, i,
                nodeID_PtrTag, MPI_COMM_WORLD,
                &(requestListRecv[requestIndex + 1]));
      requestIndex += 2;
    }
  }
  // Check if sends tally with element recv count
  assert(requestIndex == 2 * recvCount);

  free(elemID_recvCount);
  free(elemID_recvCountCum);
  free(nodePtr);
  free(nodeCount);
  free(nodeList);
  free(nodeID_countSend);

  // Receive the node IDs
  requestIndex = 0;
  for (int i = 0; i < world_size; ++i) {
    if (ghostCountProcess[i]) {
      int location = nodeID_countRecvCum[i];
      MPI_Irecv(&(nodeListRecv[location]), nodeID_countRecv[i], MPI_INT, i,
                nodeIDTag, MPI_COMM_WORLD, &(requestListSend[requestIndex]));
      MPI_Irecv(&(nodePtrRecv[uniqueElementCountCum[i] + 1]),
                ghostCountProcess[i], MPI_INT, i, nodeID_PtrTag, MPI_COMM_WORLD,
                &(requestListSend[requestIndex + 1]));
      requestIndex += 2;
    }
  }
  assert(requestIndex == 2 * sendCount);

  // Wait for all sends and recvs to be over
  for (int i = 0; i < 2 * recvCount; ++i) {
    MPI_Wait(&(requestListRecv[i]), &status);
  }
  for (int i = 0; i < 2 * sendCount; ++i) {
    MPI_Wait(&(requestListSend[i]), &status);
  }
  free(requestListSend);
  free(requestListRecv);

  // Convert nodePtrRecv to cumulative array
  for (int i = 0; i < nodePtrSize; ++i) {
    nodePtrRecv[i + 1] += nodePtrRecv[i];
  }
#ifdef DEBUG
  if (debug && 1 == 0) {
    for (int i = 0; i < world_size; ++i) {
      if (ghostCountProcess[i]) {
        printf("DEBUG(%d) : Recv from %d\t : ", world_rank, i);
        for (int j = 0; j < ghostCountProcess[i]; ++j) {
          printf("%d\t", nodePtrRecv[uniqueElementCountCum[i] + 1 + j]);
        }
        printf("\n");
      }
    }
  }
  if (debug && 1 == 0) {
    printf("\nDEBUG(%d) : Printing available ghost data\n", world_rank);
    for (int i = 0; i < world_size; ++i) {
      printf("\nElement ID from process %d\n", i);
      for (int j = ghostCountProcessCum[i], l = 0;
           j < (ghostCountProcessCum[i] + ghostCountProcess[i]); ++j, ++l) {
        printf("%d\t : ", elemID_request[j]);
        int loc = uniqueElementCountCum[i] + l;
        for (int k = nodePtrRecv[loc]; k < nodePtrRecv[loc + 1]; ++k) {
          printf("%d\t", nodeListRecv[k]);
        }
        printf("\n");
      }
    }
    printf("\n");
  }
#endif //DEBUG
  free(ghostCountProcessCum);
  free(elemID_request);
  free(uniqueElementCountCum);
  free(nodePtrRecv);

  // Find common nodes between boundary elements and ghost elements
  // Create cumulative array of node list for each process
  int *nodeListCountProcessCum = (int *)malloc((world_size + 1) * sizeof(int));
  int *nodeListCurrentProcess = (int *)calloc(world_size, sizeof(int));
  nodeListCountProcessCum[0] = 0;
  for (int i = 0; i < world_size; ++i) {
    nodeListCountProcessCum[i + 1] =
        nodeListCountProcessCum[i] + nodeListCountProcess[i];
  }
  int totalLocalNodes = nodeListCountProcessCum[world_size];
  int *localBoundaryNodeList = (int *)malloc(totalLocalNodes * sizeof(int));

  // Find total number of processes to communicate
  sendProcessCount = 0;
  for (int i = 0; i < world_size; ++i) {
    if (ghostCountProcess[i]) {
      sendProcessCount += 1;
    }
  }
  // Create list of all processes to communicate
  sendProcessID = (int *)malloc(sendProcessCount * sizeof(int));
  int loc = 0;
  for (int i = 0; i < world_size; ++i) {
    if (ghostCountProcess[i]) {
      sendProcessID[loc] = i;
      loc = loc + 1;
    }
  }
  // Check if all prcoess ID are assigned
  assert(loc == sendProcessCount);

  FILE_LOGArrayInt(DEBUGLOGIGNORE, sendProcessID, sendProcessCount, "List of process sharing boundary");

  for (int j = 0; j < nelements; ++j) {
    if (boundaryElement[j]) {
      int startLoc = eptr[j];
      int localElementNodeSize = eptr[j + 1] - eptr[j];
      int end = xadj[j + 1];
      for (int i = xadj[j]; i < end; ++i) {
        int processID = processID_ghost[i];
        if (processID != world_rank) {
          int location = nodeListCountProcessCum[processID] +
                         nodeListCurrentProcess[processID];
          memcpy(&(localBoundaryNodeList[location]), &(connectivity[startLoc]),
                 sizeof(int) * localElementNodeSize);
          nodeListCurrentProcess[processID] += localElementNodeSize;
        }
      }
    }
  }
  // Check if node list is fully populated
  for (int i = 0; i < world_size; ++i) {
    assert(nodeListCurrentProcess[i] == nodeListCountProcess[i]);
  }
  free(processID_ghost);
  free(boundaryElement);

#ifdef DEBUG
  if (debug && 1 == 0) {
    printf("DEBUG (%d) : Shared node list before unique \n", world_rank);
    for (int i = 0; i < world_size; ++i) {
      printf("On process %d\n", i);
      for (int j = nodeListCountProcessCum[i];
           j < nodeListCountProcessCum[i + 1]; ++j) {
        printf("%d\t", localBoundaryNodeList[j]);
      }
      printf("\n");
    }
    printf("\n");
  }
#endif //DEBUG
  // Make all local node list unique
  // Use variable nodeListCurrentProcess to store unique size
  for (int i = 0; i < world_size; ++i) {
    if (ghostCountProcess[i]) {
      int location = nodeListCountProcessCum[i];
      int size = nodeListCountProcess[i];
      qsort(&(localBoundaryNodeList[location]), size, sizeof(int), compare);
      nodeListCurrentProcess[i] =
          unique(&(localBoundaryNodeList[location]), size);
    }
  }
  free(nodeListCountProcess);

#ifdef DEBUG
  if (debug && 1 == 0) {
    printf("DEBUG (%d) : Shared node list after unique \n", world_rank);
    for (int i = 0; i < world_size; ++i) {
      printf("On process %d\n", i);
      int location = nodeListCountProcessCum[i];
      int size = nodeListCurrentProcess[i];
      for (int j = 0; j < size; ++j) {
        printf("%d\t", localBoundaryNodeList[location + j]);
      }
      printf("\n");
    }
    printf("\n");
  }
  // Make node list received from processes unique
  if (debug && 1 == 0) {
    printf("DEBUG (%d) : Received node list before unique \n", world_rank);
    for (int i = 0; i < world_size; ++i) {
      printf("On process %d\n", i);
      for (int j = nodeID_countRecvCum[i]; j < nodeID_countRecvCum[i + 1];
           ++j) {
        printf("%d\t", nodeListRecv[j]);
      }
      printf("\n");
    }
    printf("\n");
  }
#endif //DEBUG
  // Make all local node list unique
  // Use variable nodeListNeighbourProcess to store unique size
  int *nodeListNeighbourProcessCount = (int *)calloc(world_size, sizeof(int));
  for (int i = 0; i < world_size; ++i) {
    if (ghostCountProcess[i]) {
      int location = nodeID_countRecvCum[i];
      int size = nodeID_countRecv[i];
      qsort(&(nodeListRecv[location]), size, sizeof(int), compare);
      nodeListNeighbourProcessCount[i] =
          unique(&(nodeListRecv[location]), size);
    }
  }
  free(nodeID_countRecv);
#ifdef DEBUG
  if (debug && 1 == 0) {
    printf("DEBUG (%d) : Received node list after unique \n", world_rank);
    for (int i = 0; i < world_size; ++i) {
      printf("On process %d\n", i);
      int location = nodeID_countRecvCum[i];
      int size = nodeListNeighbourProcessCount[i];
      for (int j = 0; j < size; ++j) {
        printf("%d\t", nodeListRecv[location + j]);
      }
      printf("\n");
    }
    printf("\n");
  }
#endif //DEBUG

  // Create intersection of local node list and node list from neighbour process
  sendNeighbourCount = (int *)malloc(sendProcessCount * sizeof(int));
  int **intersectionStorage = (int **)malloc(sendProcessCount * sizeof(int *));
  int countLocation = 0;
  for (int i = 0; i < world_size; ++i) {
    if (ghostCountProcess[i]) {
      int location1 = nodeListCountProcessCum[i];
      int size1 = nodeListCurrentProcess[i];
      int location2 = nodeID_countRecvCum[i];
      int size2 = nodeListNeighbourProcessCount[i];
      sendNeighbourCount[countLocation] = getIntersection(
          &(localBoundaryNodeList[location1]), &(nodeListRecv[location2]),
          size1, size2, &(intersectionStorage[countLocation]));
      countLocation += 1;
    }
  }
  assert(countLocation == sendProcessCount);

  free(ghostCountProcess);
  free(nodeListCountProcessCum);
  free(nodeID_countRecvCum);
  free(nodeListRecv);
  free(nodeListCurrentProcess);
  free(localBoundaryNodeList);
  free(nodeListNeighbourProcessCount);

  // Store the list of intersections
  // TODO(Anil) Multiply by ndim of count and count cum
  sendNeighbourCountCum = (int *)malloc((sendProcessCount + 1) * sizeof(int));
  sendNeighbourCountCum[0] = 0;
  for (int i = 0; i < sendProcessCount; ++i) {
    sendNeighbourCountCum[i + 1] =
        sendNeighbourCountCum[i] + sendNeighbourCount[i];
  }
  int totalUniqueNodes = sendNeighbourCountCum[sendProcessCount];
  sendNodeIndex = (int *)malloc(totalUniqueNodes * sizeof(int));
  for (int i = 0; i < sendProcessCount; ++i) {
    int location = sendNeighbourCountCum[i];
    int size = sendNeighbourCount[i];
    memcpy(&(sendNodeIndex[location]), intersectionStorage[i],
           sizeof(int) * size);
  }

  for (int i = 0; i < sendProcessCount; ++i) {
    free(intersectionStorage[i]);
  }
  free(intersectionStorage);

#ifdef DEBUG
  if (debug && 1 == 0) {
    printf("DEBUG(%d) : Node list to share \n", world_rank);
    for (int i = 0; i < sendProcessCount; ++i) {
      printf("With process %d\n", sendProcessID[i]);
      for (int j = sendNeighbourCountCum[i]; j < sendNeighbourCountCum[i + 1];
           ++j) {
        printf("%d\t", sendNodeIndex[j]);
      }
      printf("\n");
    }
    printf("\n");
  }
#endif //DEBUG

  sendNodeDisplacement = (double *)malloc(ndim*totalUniqueNodes * sizeof(double));
  recvNodeDisplacement = (double *)malloc(ndim*totalUniqueNodes * sizeof(double));
}
