#include "FemTech.h"

void updateInternalForceNeighbour(void);
void updateHGForceNeighbour(void);

void GetForce_3D() {
  // TODO(Anil) special treatment for first time step
  // Below algorithm works for n > 0
  // Following Belytschko
  // Set force_n to zero
  // TODO : Call function to update fe in cases where external traction is
  // required or call updateExternalForce BC routine before calling get 
  memcpy(f_net, fe, nDOF * sizeof(double));
  memset(fi, 0, nDOF * sizeof(double));
  memset(f_hg, 0, nDOF * sizeof(double));

  // Loop over elements and Gauss points
  for (int i = 0; i < nelements; i++) {
    int nNodesL = nShapeFunctions[i];
    // number of shape functions * ndim
    // TODO : nShapeFunctions not always equal to nNodes
    // TODO : remove the following calloc to global variables
    double *fintLocal = (double*)calloc(nNodesL*ndim, sizeof(double));
    double *f_hgLocal = (double*)calloc(nNodesL*ndim, sizeof(double));
		// force calculaton for hexes, tets and quads
		for(int j=0; j<GaussPoints[i]; j++) {
			// truss elements are unique b/c the 3D formulation is not like solid
			// elmeents like quads or hexes, or tets.
	    if(strcmp(ElementType[i], "T3D2") == 0){
				TrussStressForceUpdate(i,j,fintLocal);
				//CalculateDeformationGradient(i, j);
        // TerminateFemTech(1);
			} else {
	      // // calculate F^n
				// CalculateDeformationGradient(i, j);
	      // // Calculate Determinant of F
				// DeterminateF(i, j);
	      // // Calculate sigma^n
				// StressUpdate(i, j);
			  // InternalForceUpdate(i, j, fintLocal);

        /* Updated Lagrangian approach */
        // Calculate quantites required for Stress computation
        // Calculate F_Xi^n
        CalculateF_XiAndInverse(i, j);
        // Calculate sigma^n
	StressUpdate(i, j);
	InternalForceUpdateUL(i, j, fintLocal);
	if(strcmp(ElementType[i], "C3D8R") == 0){
          ComputeHG(i, fintLocal, f_hgLocal);
        }
			} // else
		} //loop on gauss points
    // Move Local internal for to global force
    for (int k = 0; k < nNodesL; ++k) {
      int dIndex = connectivity[eptr[i] + k];
      for (int l = 0; l < ndim; ++l) {
        fi[dIndex * ndim + l] += fintLocal[k * ndim + l];
	f_hg[dIndex * ndim + l] += f_hgLocal[k * ndim + l];
      }
    }
    free(fintLocal);
    free(f_hgLocal);
  } // loop on i, nelements
  updateInternalForceNeighbour();
  // Update net force with internal force
  for (int i = 0; i < nDOF; ++i) {
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

  FILE_LOGArray(DEBUGLOGIGNORE, fi, nDOF, "Lumped Internal Force After Exchange");
}
void updateHGForceNeighbour(void) {
  // Update array to send
  int totalNodeToSend = sendNeighbourCountCum[sendProcessCount];
  for (int i = 0; i < totalNodeToSend; ++i) {
    int nodeIndex = ndim * sendNodeIndex[i];
    memcpy(&(sendNodeDisplacement[ndim * i]), &(f_hg[nodeIndex]),
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
      f_hg[nodeIndex + l] += recvNodeDisplacement[ndim * i + l];
    }
  }
  free(requestListSend);
  free(requestListRecv);

  FILE_LOGArray(DEBUGLOGIGNORE, f_hg, nDOF, "Lumped Internal Force After Exchange");
}
