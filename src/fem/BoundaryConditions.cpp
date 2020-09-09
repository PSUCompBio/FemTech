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

  FILE_LOGMatrix_SINGLE(DEBUGLOG, stiffness, modifiedMatSize, modifiedMatSize, \
      "Printing Full Stiffness Matrix");
  int offset = nDOF*modifiedMatSize;
  FILE_LOGMatrix_SINGLE(DEBUGLOG, stiffness+offset, nSpecifiedDispBC, modifiedMatSize, \
      "Printing Auxillary Stiffness Matrix");
  FILE_LOGArraySingle(DEBUGLOG, rhs, modifiedMatSize, "Printing RHS");
  FILE_LOGArraySingle(DEBUGLOG, rhs+modifiedMatSize, nSpecifiedDispBC, "Printing Aux RHS");
}
