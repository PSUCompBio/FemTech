#include "FemTech.h"

double *rhs;
int nSpecifiedDispBC;

void ApplySteadyBoundaryConditions(void) {
  rhs = (double*)calloc(nDOF, sizeof(double));
  nSpecifiedDispBC = 0;
  for (int i = 0; i < nDOF; ++i) {
    if(boundary[i] != 0) {
      nSpecifiedDispBC += 1;
    }
  }
  // Eliminate columns
  if (nSpecifiedDispBC == 0) {
   return;
  }
  int writeColumn = 0;
  // Eliminate columns first
  for (int i = 0; i < nDOF; ++i) {
    if(boundary[i] == 0) {
      // Retain column
      if (writeColumn == i) {
        writeColumn = writeColumn + 1;
        continue;
      } else {
        double *writeLocation = &(stiffness[writeColumn*nDOF]);
        double *readLocation = &(stiffness[i*nDOF]);
        memcpy(writeLocation, readLocation, sizeof(double)*nDOF);
        writeColumn += 1;
      }
    } else {
      // Modify RHS
      double displ = displacements[i];
      for (int j = 0; j < nDOF; ++j) {
        rhs[j] -= displ*stiffness[i*nDOF+j];
      }
    }
  }
  // eliminate rows
  // Eliminated rows are stored in latter part of stiffness matrix so as to
  // computed tractions at nodes with prescribed displacement. Space is present
  // in the stiffness matrix from the eliminated columns. 
  int modifiedMatSize = nDOF-nSpecifiedDispBC;
  int rIndex = 0, wIndex = 0; // read and write index to eliminate row
  int wAuxIndex = nDOF*modifiedMatSize; //write index for eliminated rows
  for (int i = 0; i < modifiedMatSize; ++i) {
    for (int j = 0; j < nDOF; ++j) {
      if(boundary[j] == 0) {
        // Retain row
        if (rIndex == wIndex) {
          wIndex = wIndex+1;
          rIndex = rIndex+1;
        } else {
          stiffness[wIndex] = stiffness[rIndex];
          wIndex = wIndex+1;
          rIndex = rIndex+1;
        }
      } else {
        // Eliminate row and move to auxillary part
        stiffness[wAuxIndex] = stiffness[rIndex];
        wAuxIndex = wAuxIndex+1;
        rIndex = rIndex+1;
      }
    }
  }
  // Modify RHS
  rIndex = wIndex = 0; // read and write index to eliminate row
  wAuxIndex = 0; //write index for eliminated rows
  double *rhsAux = (double*)malloc(sizeof(double)*nSpecifiedDispBC);
  for (int j = 0; j < nDOF; ++j) {
    if(boundary[j] == 0) {
      // Retain row
      if (j == wIndex) {
        wIndex = wIndex+1;
      } else {
        rhs[wIndex] = rhs[j];
        wIndex = wIndex+1;
      }
    } else {
      // Eliminate row and move to auxillary part
      rhsAux[wAuxIndex] = rhs[j];
      wAuxIndex = wAuxIndex+1;
    }
  }
  // Copy aux rhs to latter part of rhs
  memcpy(&(rhs[modifiedMatSize]), rhsAux, sizeof(double)*nSpecifiedDispBC);
  free(rhsAux);

#ifdef DEBUG
  if (debug) {
    printf("DEBUG : Printing Full Stiffness Matrix\n");
    for (int j = 0; j < modifiedMatSize; ++j) {
      for (int k = 0; k < modifiedMatSize; ++k) {
        printf("%12.4f", stiffness[j+k*modifiedMatSize]);
      }
      printf("\n");
    }
    printf("DEBUG : Printing Auxillary Stiffness Matrix\n");
    int offset = nDOF*modifiedMatSize;
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
    for (int j = modifiedMatSize; j < nDOF; ++j) {
      printf("%12.4f\n", rhs[j]);
    }
  }
#endif //DEBUG
}
