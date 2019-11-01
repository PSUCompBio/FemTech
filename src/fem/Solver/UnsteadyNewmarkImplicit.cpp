#include "FemTech.h"
#include "blas.h"

// double *rhs;
// int nSpecifiedDispBC;

void ModifyMassAndStiffnessMatrix();
void ComputeUnsteadyRHS(const int n, const int nMax, double* displacementFinal, \
    const double a, const double b, const double c);

void SolveUnsteadyNewmarkImplicit(double beta, double gamma, double dt, \
    double timeFinal, char* name) {
  // Make a copy of displacements to convert it into unsteady load increasing
  // incrementally from 0 to timeFinal

  const int nDOF = nnodes*ndim;
  rhs = (double*)calloc(nDOF, sizeof(double));
  double *displacementFinal = (double*)malloc(nDOF*sizeof(double));
  memcpy(displacementFinal, displacements, sizeof(double)*nDOF);
  double *accelerationsOld = (double*)malloc(nDOF*sizeof(double));
  double *displacementsOld = (double*)malloc(nDOF*sizeof(double));

  int nMax = (int)((timeFinal+1e-7)/dt);

  // Precompute multipliers
  const double a = 1.0/(beta*dt*dt);
  const double b = 1.0/(beta*dt);
  const double c = 0.5/beta-1.0;
  const double d = dt*(1.0-gamma);
  const double e = dt*gamma;
  // Create the modified stiffness matrix
  for (int i = 0; i < nDOF*nDOF; ++i) {
    stiffness[i] += a*mass[i];
  }

  // Eliminate rows and colums of modiffied stiffness and mass matrix
  ModifyMassAndStiffnessMatrix();
  int matSize = nDOF-nSpecifiedDispBC;

  // LU decompose modified stiffness matrix
  int *pivot = (int*)malloc(matSize*sizeof(int));
  int info;
  int oneI = 1;
  dgetrf_(&matSize, &matSize, stiffness, &matSize, pivot, &info); 
  if (info) {
    printf("LU Decomposition failed with info code %d\n", info);
  }

  for (int n = 1; n < nMax+1; ++n) {
    Time = n*dt;
    printf("\nTime : %.4f\n", Time);
    // Compute the RHS for the current time step
    memcpy(accelerationsOld, accelerations, nDOF*sizeof(double));
    memcpy(displacementsOld, displacements, nDOF*sizeof(double));
    ComputeUnsteadyRHS(n, nMax, displacementFinal, a, b, c);

    dgetrs_(chn, &matSize, &oneI, stiffness, &matSize, pivot, rhs, &matSize, &info);
    if (info) {
      printf("Matrix Solution failed during back substitution with info code %d\n", info);
    }
    // Update new accelerations and velocities
    int j = 0;
    for (int i = 0; i < nDOF; ++i) {
      if(boundary[i] == 0) {
        displacements[i] = rhs[j];
        j = j + 1;
      }    
    }
    // Update velocities and displacements
    for (int i = 0; i < nDOF; ++i) {
      accelerations[i] = a*(displacements[i]-displacementsOld[i])-b*velocities[i]-c*accelerationsOld[i];
      velocities[i] += e*accelerations[i]+d*accelerationsOld[i];
    }
#ifdef DEBUG
    if (debug) {
      printf("DEBUG : Printing Displacement Solution\n");
      for (int i = 0; i < nnodes; ++i) {
        for (int j = 0; j < ndim; ++j) {
          printf("%12.4f", displacements[i*ndim+j]);
        }
        printf("\n");
      }
    }
#endif //DEBUG
    // Write the solution to file
    WriteVTU(name, n, Time);
  }
  // Free all dynamically allocated ememory
  free(displacementFinal);
  free(displacementsOld);
  free(accelerationsOld);
  free(pivot);
}

void ModifyMassAndStiffnessMatrix() {
  const int matSize = nnodes*ndim;
  nSpecifiedDispBC = 0;
  for (int i = 0; i < matSize; ++i) {
    if(boundary[i] != 0) {
      nSpecifiedDispBC += 1;
    }
  }
#ifdef DEBUG
  printf("Debug : %d\n", nSpecifiedDispBC);
#endif //DEBUG
  // Eliminate columns
  if (nSpecifiedDispBC == 0) {
   return;
  }
  // Allocate space to save eliminated columns, to use in computing modified rhs
  double *stiffnessEliminated = (double*)malloc(matSize*nSpecifiedDispBC*sizeof(double));
  int writeColumn = 0;
  int writeAuxColumn = 0;
  // Eliminate columns first
  for (int i = 0; i < matSize; ++i) {
    if(boundary[i] == 0) {
      // Retain column
      if (writeColumn == i) {
        writeColumn = writeColumn + 1;
        continue;
      } else {
        double *writeLocation = &(stiffness[writeColumn*matSize]);
        double *readLocation = &(stiffness[i*matSize]);
        memcpy(writeLocation, readLocation, sizeof(double)*matSize);
        writeColumn += 1;
      }
    } else {
      double *writeAuxLocation = &(stiffnessEliminated[writeAuxColumn*matSize]);
      double *readLocation = &(stiffness[i*matSize]);
      memcpy(writeAuxLocation, readLocation, sizeof(double)*matSize);
      writeAuxColumn += 1;
    }
  }
  // Move eliminated columns to end of matrix
  const int modifiedMatSize = matSize-nSpecifiedDispBC;
  memcpy(&stiffness[matSize*modifiedMatSize], stiffnessEliminated, \
        matSize*nSpecifiedDispBC*sizeof(double));
  // eliminate rows
  // Eliminated rows are stored in middle part of stiffness matrix so as to
  // computed tractions at nodes with prescribed displacement.  
  int rIndex = 0, wIndex = 0; // read and write index to eliminate row
  int wAuxIndex = 0; //write index for eliminated rows
  for (int i = 0; i < modifiedMatSize; ++i) {
    for (int j = 0; j < matSize; ++j) {
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
        stiffnessEliminated[wAuxIndex] = stiffness[rIndex];
        wAuxIndex = wAuxIndex+1;
        rIndex = rIndex+1;
      }
    }
  }
  // Move eliminated row to middle of matrix
  memcpy(&stiffness[modifiedMatSize*modifiedMatSize], stiffnessEliminated, \
        modifiedMatSize*nSpecifiedDispBC*sizeof(double));

  // Mass matrix is used in RHS and to compute final traction
  // Hence eliminate olny rows and store at the end of matrix
  rIndex = 0, wIndex = 0; // read and write index to eliminate row
  wAuxIndex = 0; //write index for eliminated rows
  for (int i = 0; i < matSize; ++i) {
    for (int j = 0; j < matSize; ++j) {
      if(boundary[j] == 0) {
        // Retain row
        if (rIndex == wIndex) {
          wIndex = wIndex+1;
          rIndex = rIndex+1;
        } else {
          mass[wIndex] = mass[rIndex];
          wIndex = wIndex+1;
          rIndex = rIndex+1;
        }
      } else {
        // Eliminate row and move to auxillary part
        // stiffnessEliminated is used as a temperory variable
        stiffnessEliminated[wAuxIndex] = mass[rIndex];
        wAuxIndex = wAuxIndex+1;
        rIndex = rIndex+1;
      }
    }
  }
  // Move eliminated row to middle of mass matrix
  memcpy(&mass[matSize*modifiedMatSize], stiffnessEliminated, \
        matSize*nSpecifiedDispBC*sizeof(double));
  // Free temperory storage
  free(stiffnessEliminated);
  return;
}

void ComputeUnsteadyRHS(const int n, const int nMax, double* displacementFinal, \
    const double a, const double b, const double c) {
  int nDOF = nnodes*ndim;
  int nRHS = nDOF-nSpecifiedDispBC;
  // Set rhs to zero
  memset(rhs, 0, sizeof(double)*nDOF);

  // Include effect from values from previous time steps
  double *rhsTn = (double*)malloc(nDOF*sizeof(double));
  for (int i = 0; i < nDOF; ++i) {
    rhsTn[i] = a*displacements[i]+b*velocities[i]+c*accelerations[i];
  }
  // Compute M times rhsTn
  double *rhs2 = (double*)malloc(nRHS*sizeof(double));
  dgemv_(chn, &nRHS, &nDOF, &one, mass, &nRHS, rhsTn, &oneI, &zero, rhs2, &oneI);
  // Set displacement at specified dofs from bc
  for (int i = 0; i < nDOF; ++i) {
    if(boundary[i] == 1) {
      displacements[i] = ((double)n/(double)nMax)*displacementFinal[i];
    }
  }
  // Include effect of eliminated columns
  int k = 0;
  double *stiffnessEliminated = &(stiffness[nDOF*nRHS]);
  for (int i = 0; i < nDOF; ++i) {
    if(boundary[i] == 1) {
      double displ = displacements[i];
      for (int j = 0; j < nDOF; ++j) {
        rhs[j] -= displ*stiffnessEliminated[k*nDOF+j];
      }
      k += 1;
    }
  }
  // Eliminated known DOF and store to auxRHS for traction computation
  int rIndex = 0, wIndex = 0; // read and write index to eliminate row
  int wAuxIndex = 0; //write index for eliminated rows
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
  memcpy(&(rhs[nRHS]), rhsAux, sizeof(double)*nSpecifiedDispBC);
  free(rhsAux);

  // Add rhs and contribution from previous time step
  for (int i = 0; i < nRHS; ++i) {
    rhs[i] += rhs2[i];
  }

  free(rhsTn);
  free(rhs2);
  return;
}
