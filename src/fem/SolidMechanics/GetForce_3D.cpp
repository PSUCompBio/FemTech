#include "FemTech.h"

void updateInternalForceNeighbour(void);

void GetForce_3D() {
  // TODO(Anil) special treatment for first time step
  // Below algorithm works for n > 0
  const int nDOF = nnodes * ndim;
  // Store internal force from previous step to compute energy
  memcpy(fi_prev, fi, nDOF * sizeof(double));
  // Following Belytschko
  // Set force_n to zero
  memset(f_net, 0, nDOF * sizeof(double));
  memset(fi, 0, nDOF * sizeof(double));
  // Loop over elements and Gauss points
  for (int i = 0; i < nelements; i++) {
    int nNodes = nShapeFunctions[i];
    // number of shape functions * ndim
    double *fintLocal = (double *)calloc(nNodes * ndim, sizeof(double));
    for (int j = 0; j < GaussPoints[i]; j++) {
      // calculate F^n
      CalculateDeformationGradient(i, j);
      // Calculate Determinant of F
      DeterminateF(i, j);
      // Calculate sigma^n
      StressUpdate(i, j);
      InternalForceUpdate(i, j, fintLocal);
    } // loop on gauss points
    // Move Local internal for to global force
    for (int k = 0; k < nNodes; ++k) {
      int dIndex = connectivity[eptr[i] + k];
      for (int l = 0; l < ndim; ++l) {
        fi[dIndex * ndim + l] += fintLocal[k * ndim + l];
        if (l ==2)
        printf("e : %d, node : %d, dim : %d, nID, %d, Value : %12.6f\n", i, k, l, dIndex, fintLocal[k*ndim+l]);
      }
    }
    free(fintLocal);
  } // loop on i, nelements
  updateInternalForceNeighbour();
  // Update net force with internal force
  for (int i = 0; i < nnodes * ndim; ++i) {
    f_net[i] -= fi[i];
  }
  return;
}
void updateInternalForceNeighbour(void) {
  // Update array to send
  int totalNodeToSend = sendNeighbourCountCum[sendProcessCount];
  for (int i = 0; i < totalNodeToSend; ++i) {
    int nodeIndex = ndim * sendNodeIndex[i];
    memcpy(&(sendNodeDisplacement[ndim * i]), &(fi[nodeIndex]),
           sizeof(double) * ndim);
  }
  // Create send requests
  int forceTag = 2169;
  MPI_Request *requestListSend =
      (MPI_Request *)malloc(sizeof(MPI_Request) * sendProcessCount);
  for (int i = 0; i < sendProcessCount; ++i) {
    int process = sendProcessID[i];
    int location = sendNeighbourCountCum[i] * ndim;
    int size = ndim * sendNeighbourCount[i];
    MPI_Isend(&(sendNodeDisplacement[location]), size, MPI_DOUBLE, process,
              forceTag, MPI_COMM_WORLD, &(requestListSend[i]));
  }
  // Create recv requests
  MPI_Request *requestListRecv =
      (MPI_Request *)malloc(sizeof(MPI_Request) * sendProcessCount);
  for (int i = 0; i < sendProcessCount; ++i) {
    int process = sendProcessID[i];
    int location = sendNeighbourCountCum[i] * ndim;
    int size = ndim * sendNeighbourCount[i];
    MPI_Irecv(&(recvNodeDisplacement[location]), size, MPI_DOUBLE, process,
              forceTag, MPI_COMM_WORLD, &(requestListRecv[i]));
  }
  // Wait for completion of all requests
  MPI_Status status;
  for (int i = 0; i < sendProcessCount; ++i) {
    MPI_Wait(&(requestListSend[i]), &status);
  }
  for (int i = 0; i < sendProcessCount; ++i) {
    MPI_Wait(&(requestListRecv[i]), &status);
  }
  // Update Internal Force values
  for (int i = 0; i < totalNodeToSend; ++i) {
    int nodeIndex = ndim * sendNodeIndex[i];
    for (int l = 0; l < ndim; ++l) {
      fi[nodeIndex + l] += recvNodeDisplacement[ndim * i + l];
    }
  }
  free(requestListSend);
  free(requestListRecv);
  if (debug && 1 == 0) {
    const int nDOF = nnodes * ndim;
    printf("Lumped Internal Force After Exchange\n");
    for (int j = 0; j < nDOF; ++j) {
      printf("%d  %12.6f\n", j, fi[j]);
    }
  }
}
