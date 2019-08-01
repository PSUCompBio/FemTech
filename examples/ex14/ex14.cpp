#include "FemTech.h"
#include "blas.h"

#include <assert.h>

/*Delare Functions*/
void CustomPlot();
void InitCustomPlot();
void InitBoundaryCondition(double *aMax);
void ApplyAccBoundaryConditions();

/* Global Variables/Parameters  - could be moved to parameters.h file?  */
double Time;
int nStep;
int nSteps;
int nPlotSteps = 100;
bool ImplicitStatic = false;
bool ImplicitDynamic = false;
bool ExplicitDynamic = true;
double ExplicitTimeStepReduction = 0.8;
double FailureTimeStep = 1e-11;

/* Global variables used only in this file */
const double sphereRadius = 0.08;
int nodeIDtoPlot;
bool rankForCustomPlot;
/* Global variables for bc */
int* boundaryID = NULL;
int boundarySize;
double aLin[3], bLin[3];
double aAng[3], bAng[3];
double peakTime, maxTime;

int main(int argc, char **argv) {
  double accMax[3] = {1.0, 0.0, 0.0};
  peakTime = 0.001;
  maxTime = 0.005;
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);
  // Get the number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  // Get the rank of the process
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  if (ReadInputFile(argv[1])) {
    PartitionMesh();
  }

  AllocateArrays();
  ReadMaterials();
  InitCustomPlot();
  InitBoundaryCondition(accMax);

  /* Write inital, undeformed configuration*/
  Time = 0.0;
  nStep = 0;
  WriteVTU(argv[1], nStep, Time);
  CustomPlot();

  // Dynamic Explcit solution using....
  double dt = 0.0;
  double tMax = 0.01;   // max simulation time in seconds

  double Time = 0.0;
  int time_step_counter = 0;
  int plot_counter = 0;
  const int nDOF = nnodes * ndim;
  /** Central Difference Method - Beta and Gamma */
  // double beta = 0;
  // double gamma = 0.5;

  ShapeFunctions();
  /*  Step-1: Calculate the mass matrix similar to that of belytschko. */
  AssembleLumpedMass();

  // Used if initial velocity and acceleration BC is to be set.

  ApplyAccBoundaryConditions();
  /* Step-2: getforce step from Belytschko */
  GetForce(); // Calculating the force term.

  /* Obtain dt, according to Belytschko dt is calculated at end of getForce */
  dt = ExplicitTimeStepReduction * StableTimeStep();

  /* Step-3: Calculate accelerations */
  CalculateAccelerations();

  nSteps = (int)(tMax / dt);
  int nsteps_plot = (int)(nSteps / nPlotSteps);

  if (world_rank == 0) {
    printf("inital dt = %3.3e, nSteps = %d, nsteps_plot = %d\n", dt, nSteps,
           nsteps_plot);
  }

  time_step_counter = time_step_counter + 1;
  double t_n = 0.0;

  if (world_rank == 0) {
    printf(
        "------------------------------- Loop ----------------------------\n");
    printf("Time : %f, tmax : %f\n", Time, tMax);
  }

  /* Step-4: Time loop starts....*/
  while (Time < tMax) {
    double t_n = Time;
    double t_np1 = Time + dt;
    Time = t_np1; /*Update the time by adding full time step */
    if (world_rank == 0) {
      printf("Time : %15.6e, dt=%15.6e, tmax : %15.6e\n", Time, dt, tMax);
    }
    double dt_nphalf = dt;                 // equ 6.2.1
    double t_nphalf = 0.5 * (t_np1 + t_n); // equ 6.2.1

    /* Step 5 from Belytschko Box 6.1 - Update velocity */
    for (int i = 0; i < nDOF; i++) {
      if (boundary[i]) {
        velocities_half[i] = velocities[i];
      } else {
        velocities_half[i] =
            velocities[i] + (t_nphalf - t_n) * accelerations[i];
      }
    }

    // Store old displacements and accelerations for energy computation
    memcpy(displacements_prev, displacements, nDOF * sizeof(double));
    memcpy(accelerations_prev, accelerations, nDOF * sizeof(double));
    // Store internal external force from previous step to compute energy
    memcpy(fi_prev, fi, nDOF * sizeof(double));
    memcpy(fe_prev, fe, nDOF * sizeof(double));

    for (int i = 0; i < nDOF; i++) {
      if (!boundary[i]) {
        displacements[i] = displacements[i] + dt_nphalf * velocities_half[i];
      }
    }
    /* Step 6 Enforce boundary Conditions */
    ApplyAccBoundaryConditions();

    /* Step - 8 from Belytschko Box 6.1 - Calculate net nodal force*/
    GetForce(); // Calculating the force term.
    /* Step - 9 from Belytschko Box 6.1 - Calculate Accelerations */
    CalculateAccelerations(); // Calculating the new accelerations from total
                              // nodal forces.

    /** Step- 10 - Second Partial Update of Nodal Velocities */
    for (int i = 0; i < nDOF; i++) {
      if (!boundary[i]) {
        velocities[i] =
            velocities_half[i] + (t_np1 - t_nphalf) * accelerations[i];
      }
    }

    /** Step - 11 Checking* Energy Balance */
    CheckEnergy(Time);

    if (time_step_counter % nsteps_plot == 0) {
      plot_counter = plot_counter + 1;
      // printf("Plot %d/%d: dt=%3.2e s, Time=%3.2e s, Tmax=%3.2e s on rank :
      // %d\n", 	plot_counter,nPlotSteps,dt,Time,tMax, world_rank);
      for (int i = 0; i < nelements; i++) {
        for (int l = 0; l < ndim * ndim; l++) {
          Favg[i * ndim * ndim + l] = 0.0;
        } // initializing avg def gradient to zero for each time step
        for (int j = 0; j < GaussPoints[i]; j++) {
          SumOfDeformationGradient(i, j);
        } // calculating sum of deformation gradient for all gauss points
        for (int k = 0; k < ndim * ndim; k++) {
          Favg[i * ndim * ndim + k] =
              Favg[i * ndim * ndim + k] / GaussPoints[i];
        } // dividing by number of gauss points to get average deformation
          // gradient
        CalculateStrain(i);
      } // calculating avergae strain for every element
      printf("------Plot %d: WriteVTU by rank : %d\n", plot_counter,
             world_rank);
      WriteVTU(argv[1], plot_counter, Time);
      CustomPlot();

#ifdef DEBUG
      if (debug) {
        printf("DEBUG : Printing Displacement Solution\n");
        for (int i = 0; i < nnodes; ++i) {
          for (int j = 0; j < ndim; ++j) {
            printf("%15.6E", displacements[i * ndim + j]);
          }
          printf("\n");
        }
      }
#endif // DEBUG
    }
    time_step_counter = time_step_counter + 1;
    dt = ExplicitTimeStepReduction * StableTimeStep();
    // Barrier not a must
    MPI_Barrier(MPI_COMM_WORLD);

    nStep = plot_counter;
    // Write out the last time step
    CustomPlot();
  } // end explcit while loop
#ifdef DEBUG
  if (debug) {
    printf("DEBUG : Printing Displacement Solution\n");
    for (int i = 0; i < nnodes; ++i) {
      for (int j = 0; j < ndim; ++j) {
        printf("%15.6E", displacements[i * ndim + j]);
      }
      printf("\n");
    }
  }
#endif // DEBUG

  /* Below are things to do at end of program */
  if (world_rank == 0) {
    WritePVD(argv[1], nStep, Time);
  }
  FreeArrays();
  // Free local boundary related arrays
  if (boundaryID) {
    free(boundaryID);
  }
  MPI_Finalize();
  return 0;
}

void ApplyAccBoundaryConditions() {
  double linAcc[3], angAcc[3];
  double linVel[3], angVel[3];
  double linDispl[3], angDispl[3];
  double cm[3];

  // Compute accelerations, velocities and displacements
  // Compute angular accelerations, angular velocities and angles
  if (Time < peakTime) {
    for (int i = 0; i < ndim; ++i) {
      linAcc[i] = aLin[i]*Time;
      linVel[i] = 0.5*aLin[i]*Time*Time;
      linDispl[i] = aLin[i]*Time*Time*Time/6.0;

      // angAcc[i] = 0.0;
      // angVel[i] = 0.0;
      // angDispl[i] = 0.0;
    }
  } else {
    if (Time < maxTime) {
      for (int i = 0; i < ndim; ++i) {
        linAcc[i] = (aLin[i]+bLin[i])*peakTime-bLin[i]*Time;
        linVel[i] = (aLin[i]+bLin[i])*(peakTime*Time-0.5*peakTime*peakTime)-\
                    0.5*bLin[i]*Time*Time;
        linDispl[i] = 0.5*(aLin[i]+bLin[i])*peakTime*(peakTime*peakTime/3.0-\
            peakTime*Time+Time*Time)-bLin[i]*Time*Time*Time/6.0;
      }
    } else {
      for (int i = 0; i < ndim; ++i) {
        linAcc[i] = 0.0;
        linVel[i] = (aLin[i]+bLin[i])*(peakTime*maxTime-0.5*peakTime*peakTime)-\
                    0.5*bLin[i]*maxTime*maxTime;
        linDispl[i] = 0.5*(aLin[i]+bLin[i])*peakTime*(peakTime*peakTime/3.0-\
            peakTime*maxTime+maxTime*maxTime)-bLin[i]*maxTime*maxTime*maxTime/6.0+\
                      linVel[i]*(Time-maxTime);
      }
    }
  }
  GetBodyCenterofMass(cm);
  for (int i = 0; i < boundarySize; i++) {
    int index = boundaryID[i]*ndim;
    for (int j = 0; j < ndim; ++j) {
      displacements[index+j] = linDispl[j];
      velocities[index+j] = linVel[j];
      // For energy computations
      accelerations[index+j] = linAcc[j];
    }
  }
  if (world_rank == 0) {
    printf("INFO(%d) : Applied acceleration = (%10.5e, %10.5e, %10.5e)\n", \
        world_rank, linAcc[0], linAcc[1], linAcc[2]);
  }
  return;
}

void InitCustomPlot() {
  double xPlot = sphereRadius;
  double yPlot = 0.00;
  double zPlot = 0.00;
  double tol = 1e-5;
  FILE *datFile;
  rankForCustomPlot = false;
  int index;
  const int x = 0;
  const int y = 1;
  const int z = 2;

  for (int i = 0; i < nnodes; ++i) {
    index = i*ndim;
    if (fabs(coordinates[index + x] - xPlot) < tol &&
        fabs(coordinates[index + y] - yPlot) < tol &&
        fabs(coordinates[index + z] - zPlot) < tol) {
      nodeIDtoPlot = i;
      rankForCustomPlot = true;
      printf("INFO(%d) : nodeID for plot : %d\n", world_rank, nodeIDtoPlot);
      break;
    }
  }
  if (rankForCustomPlot) {
    datFile = fopen("plot.dat", "w");
    fprintf(datFile, "# Results for Node %d\n", nodeIDtoPlot);
    fprintf(datFile, "# Time  DispX    DispY   DispZ\n");
    fclose(datFile);
  }
  return;
}

void CustomPlot() {
  if (rankForCustomPlot) {
    const int x = 0;
    const int y = 1;
    const int z = 2;
    const int index = nodeIDtoPlot*ndim;

    FILE *datFile;
    datFile = fopen("plot.dat", "a");
    fprintf(datFile, "%11.3e %11.3e  %11.3e  %11.3e\n", Time,
            displacements[index + x], displacements[index + y],
            displacements[index + z]);

    fclose(datFile);
  }
  return;
}

void InitBoundaryCondition(double* aMax) {
  double tol = 1e-5;
  // Find the number of nodes on the outer boundary
  // For sphere : find points at specified radius
  double sphereRad2 = sphereRadius*sphereRadius;
  int index;
  double nodeDist2;
  int boundarySize = 0;
  for (int i = 0; i < nnodes; ++i) {
    index = i*ndim;
    nodeDist2 = coordinates[index]*coordinates[index]+coordinates[index+1]*\
                coordinates[index+1]+coordinates[index+2]*coordinates[index+2];
    if (fabs(nodeDist2-sphereRad2) < tol) {
      boundarySize = boundarySize + 1;
    }
  }
  if (boundarySize) {
    boundaryID = (int*)malloc(boundarySize*sizeof(int));
    if (boundaryID == NULL) {
      printf("ERROR(%d) : Unable to alocate boundaryID\n", world_rank);
      exit(0);
    }
  }
  int bcIndex = 0;
  for (int i = 0; i < nnodes; ++i) {
    index = i*ndim;
    nodeDist2 = coordinates[index]*coordinates[index]+coordinates[index+1]*\
                coordinates[index+1]+coordinates[index+2]*coordinates[index+2];
    if (fabs(nodeDist2-sphereRad2) < tol) {
      boundaryID[bcIndex] = i;
      bcIndex = bcIndex + 1;
    }
  }
  assert(bcIndex == boundarySize);
  // Compute the constants required for acceleration computations
  for (int i = 0; i < ndim; ++i) {
    aLin[i] = aMax[i]/peakTime;
    bLin[i] = aMax[i]/(maxTime-peakTime);
  }
}
