#include "FemTech.h"

#include "blas.h"

void MassElementMatrix(double *Me, int e) {
  // number of shape functions * ndim
  int bColSize = nShapeFunctions[e]*ndim;
  // Size of local mass matrix [nShapeFunctions*ndim x nShapeFunctions*ndim]
  int mLocalSize = bColSize*bColSize;

  //TODO(Anil) remove allocations inside the loop
  double *MeGQ = (double*)calloc(mLocalSize, sizeof(double));
  // Column size of N matrix
  int nColSize = nShapeFunctions[e]*ndim;
  // Size of N matrix [ndim x nNodes*ndim]
  int Nsize = nColSize*ndim;
  int NiSize = ndim*ndim;
  int NindexStart;
  double Ni;
  // Create N Matrix of size ndim x nNodes*ndim in column major format.
  // N = [N_1, N_2, N_3, ... N_nNodes]
  // N_n = [N_i 0 0
  //        0 N_i 0
  //        0 0 N_i]
  // N is computed for each Gauss Quadrature point
  double *N = (double*)calloc(Nsize, sizeof(double));
  for (int k = 0; k < GaussPoints[e]; k++) {
    // Populate N for each Gauss Point
    for (int n = 0; n < nShapeFunctions[e]; ++n) {
      NindexStart = n*NiSize;
      Ni = shp[gptr[e]+k*GaussPoints[e]+n];

      // N_n = [N_i 0 0
      //        0 N_i 0
      //        0 0 N_i]
      // N_n stored in column major format N(1, 1) = N[0], N(2, 2) = N[4], N(3, 3) = N[8];
      N[NindexStart] = Ni;
      N[NindexStart+4] = Ni;
      N[NindexStart+8] = Ni;
    }
    // Compute elemental mass metrix for a single Gauss Quadrature point MeGQ = N^T N
    dgemm_(chy, chn, &nColSize, &nColSize, &ndim, &one, N, &ndim, \
        N, &ndim, &zero, MeGQ, &nColSize);
    // TODO(Anil) Gauss weights and det J product can be combined in shape
    // function
    int wIndex = gpPtr[e]+k;
    const double preFactor = gaussWeights[wIndex]*detJacobian[wIndex];
    // Me = \Sum_j w_j (N^T*N Det(J))_j
    // Add contribution of current Gauss point to elemental mass matrix
    for (int n = 0; n < mLocalSize; ++n) {
      Me[n] += MeGQ[n]*preFactor;
    }
  }
  // Multiply all elements of mass matrix by density
  // TODO(Anil) Use material property rho and remove variable rho
  for (int n = 0; n < mLocalSize; ++n) {
    Me[n] *= rho;
  }
  // print Me Matrix
  if (debug && 1==1) {
    printf("DEBUG : Printing Me (Mass Matrix) for Element %d\n", e);
    for (int j = 0; j < bColSize; ++j) {
      for (int k = 0; k < bColSize; ++k) {
        printf("%.4f\t", Me[j+k*bColSize]);
      }
      printf("\n");
    }
  }
  free(MeGQ);
  free(N);
  return;
}
void LumpMassMatrix(void) {
  const int massSize = nnodes*ndim;
  // Lump mass matrix by summing up along the row
  for(int i = 1; i < massSize; ++i) {
    for(int j = 0; j < massSize; ++j) {
      mass[j] += mass[j+i*massSize];
    }
  }
	if(debug && 1==0){
	  printf("Lumped Mass\n");
	  for(int j = 0; j < massSize; ++j) {
	    printf("%d  %12.6f\n", j, mass[j]);
	  }
	}
}
void updateMassMatrixNeighbour(void) {
  // Update array to send 
  int totalNodeToSend = sendNeighbourCountCum[sendProcessCount];
  for (int i = 0; i < totalNodeToSend; ++i) {
    int nodeIndex = ndim*sendNodeIndex[i];
    sendNodeDisplacement[ndim*i] = mass[nodeIndex];
    sendNodeDisplacement[ndim*i+1] = mass[nodeIndex+1];
    sendNodeDisplacement[ndim*i+2] = mass[nodeIndex+2];
  }
  // Create send requests
  int massTag = 2168;
  MPI_Request *requestListSend =
        (MPI_Request *)malloc(sizeof(MPI_Request) * sendProcessCount);
  for (int i = 0; i < sendProcessCount; ++i) {
    int process = sendProcessID[i];
    int location = sendNeighbourCountCum[i]*ndim;
    int size = ndim*sendNeighbourCount[i];
    MPI_Isend(&(sendNodeDisplacement[location]), size, MPI_DOUBLE, process, \
        massTag, MPI_COMM_WORLD, &(requestListSend[i]));
  }
  // Create recv requests
  MPI_Request *requestListRecv =
        (MPI_Request *)malloc(sizeof(MPI_Request) * sendProcessCount);
  for (int i = 0; i < sendProcessCount; ++i) {
    int process = sendProcessID[i];
    int location = sendNeighbourCountCum[i]*ndim;
    int size = ndim*sendNeighbourCount[i];
    MPI_Irecv(&(recvNodeDisplacement[location]), size, MPI_DOUBLE, process, \
        massTag, MPI_COMM_WORLD, &(requestListRecv[i]));
  }
  // Wait for completion of all requests
  MPI_Status status;
  for (int i = 0; i < sendProcessCount; ++i) {
    MPI_Wait(&(requestListSend[i]), &status);
  }
  for (int i = 0; i < sendProcessCount; ++i) {
    MPI_Wait(&(requestListRecv[i]), &status);
  }
  // Update Mass values
  for (int i = 0; i < totalNodeToSend; ++i) {
    int nodeIndex = ndim*sendNodeIndex[i];
    mass[nodeIndex] += recvNodeDisplacement[ndim*i];
    mass[nodeIndex+1] += recvNodeDisplacement[ndim*i+1];
    mass[nodeIndex+2] += recvNodeDisplacement[ndim*i+2];
  }
  free(requestListSend);
  free(requestListRecv);
	if(debug && 1==0){
    const int massSize = nnodes*ndim;
	  printf("Lumped Mass After Exchange\n");
	  for(int j = 0; j < massSize; ++j) {
	    printf("%d  %12.6f\n", j, mass[j]);
	  }
	}
}
