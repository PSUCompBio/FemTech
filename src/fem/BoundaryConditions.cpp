#include "FemTech.h"

double *rhs;

void ApplySteadyBoundaryConditions(void) {
  const int matSize = nnodes*ndim;
  rhs = (double*)calloc(matSize, sizeof(double));
  // Eliminate columns
  if (nSpecifiedDispBC == 0) {
   return;
  }
  // First column to eliminate
  int colCount = nodeDOFDispBC[0];
  for (int i = 0; i < nSpecifiedDispBC-1; ++i) {
    // Update RHS corresponding to eliminated column
    int colNumber = nodeDOFDispBC[i];
    double displ = nodeValueDispBC[i];
    for (int j = 0; j < matSize; ++j) {
      rhs[j] -= displ*stiffness[colNumber*matSize+j];
    }
    // Move columns to left eliminating current column
    for(int j = colNumber+1; j < nodeDOFDispBC[i+1]; ++j) {
      double *writeLocation = &(stiffness[colCount*matSize]);
      double *readLocation = &(stiffness[j*matSize]);
      memcpy(writeLocation, readLocation, sizeof(double)*matSize);
      colCount += 1;
    }
  }
  // Update RHS to eliminated last column
  int colNumber = nodeDOFDispBC[nSpecifiedDispBC-1];
  double displ = nodeValueDispBC[nSpecifiedDispBC-1];
  for (int j = 0; j < matSize; ++j) {
    rhs[j] -= displ*stiffness[colNumber*matSize+j];
  }
  // Overwrite from last column to eliminate to end of matrix
  for(int j = nodeDOFDispBC[nSpecifiedDispBC-1]; j < matSize; ++j) {
    double *writeLocation = &(stiffness[colCount*matSize]);
    double *readLocation = &(stiffness[j*matSize]);
    memcpy(writeLocation, readLocation, sizeof(double)*matSize);
    colCount += 1;
  }
  // eliminate rows
  // Eliminated rows are stored in latter part of stiffness matrix so as to
  // computed tractions at nodes with prescribed displacement. Space is present
  // in the stiffness matrix from the eliminated columns. 
  int modifiedMatSize = matSize-nSpecifiedDispBC;
  int rIndex = 0, wIndex = 0; // read and write index to eliminate row
  int wAuxIndex = matSize*modifiedMatSize; //write index for eliminated rows
  for (int i = 0; i < modifiedMatSize; ++i) {
    for (int j = 0; j < nodeDOFDispBC[0]; ++j) {
      stiffness[wIndex] = stiffness[rIndex];
      wIndex = wIndex+1;
      rIndex = rIndex+1;
    }
    for (int k = 0; k < nSpecifiedDispBC-1; ++k) {
      int node = nodeDOFDispBC[k];
      stiffness[wAuxIndex] = stiffness[rIndex];
      wAuxIndex = wAuxIndex+1;
      rIndex = rIndex+1;
      for (int j = node+1; j < nodeDOFDispBC[k+1]; ++j) {
        stiffness[wIndex] = stiffness[rIndex];
        wIndex = wIndex+1;
        rIndex = rIndex+1;
      }
    }
    int node = nodeDOFDispBC[nSpecifiedDispBC-1];
    stiffness[wAuxIndex] = stiffness[rIndex];
    wAuxIndex = wAuxIndex+1;
    rIndex = rIndex+1;
    for (int j = node+1; j < matSize; ++j) {
      stiffness[wIndex] = stiffness[rIndex];
      wIndex = wIndex+1;
      rIndex = rIndex+1;
    }
  }
  // Modify RHS
  rIndex = wIndex = 0; // read and write index to eliminate row
  wAuxIndex = 0; //write index for eliminated rows
  double *rhsAux = (double*)malloc(sizeof(double)*nSpecifiedDispBC);
  for (int j = 0; j < nodeDOFDispBC[0]; ++j) {
    rhs[wIndex] = rhs[rIndex];
    wIndex = wIndex+1;
    rIndex = rIndex+1;
  }
  for (int k = 0; k < nSpecifiedDispBC-1; ++k) {
    int node = nodeDOFDispBC[k];
    rhsAux[wAuxIndex] = rhs[rIndex];
    wAuxIndex = wAuxIndex+1;
    rIndex = rIndex+1;
    for (int j = node+1; j < nodeDOFDispBC[k+1]; ++j) {
      rhs[wIndex] = rhs[rIndex];
      wIndex = wIndex+1;
      rIndex = rIndex+1;
    }
  }
  int node = nodeDOFDispBC[nSpecifiedDispBC-1];
  rhsAux[wAuxIndex] = rhs[rIndex];
  wAuxIndex = wAuxIndex+1;
  rIndex = rIndex+1;
  for (int j = node+1; j < matSize; ++j) {
    rhs[wIndex] = rhs[rIndex];
    wIndex = wIndex+1;
    rIndex = rIndex+1;
  }
  // Copy aux rhs to latter part of rhs
  memcpy(&(rhs[modifiedMatSize]), rhsAux, sizeof(double)*nSpecifiedDispBC);
  free(rhsAux);
  if (debug) {
    printf("DEBUG : Printing Full Stiffness Matrix\n");
    for (int j = 0; j < modifiedMatSize; ++j) {
      for (int k = 0; k < modifiedMatSize; ++k) {
        printf("%12.4f", stiffness[j+k*modifiedMatSize]);
      }
      printf("\n");
    }
    printf("DEBUG : Printing Auxillary Stiffness Matrix\n");
    int offset = matSize*modifiedMatSize;
    int i = 0;
    for (int j = 0; j < nSpecifiedDispBC; ++j) {
      for (int k = 0; k < modifiedMatSize; ++k) {
        printf("%12.4f", stiffness[offset+j+k*nSpecifiedDispBC]);
      }
      printf("\n");
    }
    printf("DEBUG : Printing RHS\n");
    for (int j = 0; j < modifiedMatSize; ++j) {
      printf("%12.4f\n", rhs[j]);
    }
    printf("DEBUG : Printing Auxillary RHS\n");
    for (int j = modifiedMatSize; j < matSize; ++j) {
      printf("%12.4f\n", rhs[j]);
    }
  }
}
