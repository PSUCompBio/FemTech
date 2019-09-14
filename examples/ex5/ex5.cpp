#include "FemTech.h"
#include "blas.h"

#include <assert.h>

/*Delare Functions*/
void CustomPlot();
void InitCustomPlot();
void InitBoundaryCondition(double *aMax, double angMax);
void ApplyAccBoundaryConditions();

/* Global Variables/Parameters */
double Time;
int nStep;
int nSteps;
int nPlotSteps = 50;
bool ImplicitStatic = false;
bool ImplicitDynamic = false;
bool ExplicitDynamic = true;
double ExplicitTimeStepReduction = 0.8;
double FailureTimeStep = 1e-11;
static const double radToDeg = 180.0/(atan(1.0)*4.0);

/* Global variables used only in this file */
int nodeIDtoPlot;
bool rankForCustomPlot;
/* Global variables for bc */
int* boundaryID = NULL;
int boundarySize;
double aLin[3], bLin[3];
double aAng, bAng;
double peakTime, tMax;
double thetaOld = 0.0;
double linDisplOld[3];
double angNormal[3];
const int rigidPartID = 0; // part ID of elements to be made rigid

int main(int argc, char **argv) {
  double accMax[3] = {32.0*9.81, 0.0, 0.0};
  // double angAccMax = 1000.0;
  double angAccMax = 0.0;
  angNormal[0] = 0.0; angNormal[1] = 0.0; angNormal[2] = 1.0;
  peakTime = 0.020;
  tMax = 0.040;
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
  InitBoundaryCondition(accMax, angAccMax);

  /* Write inital, undeformed configuration*/
  Time = 0.0;
  nStep = 0;
  WriteVTU(argv[1], nStep, Time);
  CustomPlot();

  // Dynamic Explcit solution using....
  double dt = 0.0;

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
  double linAcc[3], angAcc;
  double linVel[3], angVel;
  double linDispl[3], angDispl;
  double cm[3];
  double position[3], rotation[3];
  double velRotation[3], omega[3];
  double alpha[3], accCorioli[3], accRotation[3], centrifugal[3];
  double rotMat[3][3];

  // Compute accelerations, velocities and displacements
  // Compute angular accelerations, angular velocities and angles
  if (Time < peakTime) {
    for (int i = 0; i < ndim; ++i) {
      linAcc[i] = aLin[i]*Time;
      linVel[i] = 0.5*aLin[i]*Time*Time;
      linDispl[i] = aLin[i]*Time*Time*Time/6.0;
    }

    angAcc = aAng*Time;
    angVel = 0.5*aAng*Time*Time;
    angDispl = aAng*Time*Time*Time/6.0;
  } else {
    if (Time < tMax) {
      for (int i = 0; i < ndim; ++i) {
        linAcc[i] = (aLin[i]+bLin[i])*peakTime-bLin[i]*Time;
        linVel[i] = (aLin[i]+bLin[i])*(peakTime*Time-0.5*peakTime*peakTime)-\
                    0.5*bLin[i]*Time*Time;
        linDispl[i] = 0.5*(aLin[i]+bLin[i])*peakTime*(peakTime*peakTime/3.0-\
            peakTime*Time+Time*Time)-bLin[i]*Time*Time*Time/6.0;
      }

      angAcc = (aAng+bAng)*peakTime-bAng*Time;
      angVel = (aAng+bAng)*(peakTime*Time-0.5*peakTime*peakTime)-\
                  0.5*bAng*Time*Time;
      angDispl = 0.5*(aAng+bAng)*peakTime*(peakTime*peakTime/3.0-\
          peakTime*Time+Time*Time)-bAng*Time*Time*Time/6.0;
    } else {
      for (int i = 0; i < ndim; ++i) {
        linAcc[i] = 0.0;
        linVel[i] = (aLin[i]+bLin[i])*(peakTime*tMax-0.5*peakTime*peakTime)-\
                    0.5*bLin[i]*tMax*tMax;
        linDispl[i] = 0.5*(aLin[i]+bLin[i])*peakTime*(peakTime*peakTime/3.0-\
            peakTime*tMax+tMax*tMax)-bLin[i]*tMax*tMax*tMax/6.0+\
                      linVel[i]*(Time-tMax);
      }
      angAcc = 0.0;
      angVel = (aAng+bAng)*(peakTime*tMax-0.5*peakTime*peakTime)-\
                  0.5*bAng*tMax*tMax;
      angDispl = 0.5*(aAng+bAng)*peakTime*(peakTime*peakTime/3.0-\
          peakTime*tMax+tMax*tMax)-bAng*tMax*tMax*tMax/6.0+\
                    angVel*(Time-tMax);
    }
  }
  double dTheta = angDispl-thetaOld;
  thetaOld = angDispl;
  GetBodyCenterofMass(cm);
  get3dRotationMatrix(angNormal, dTheta, rotMat);
  for (int j = 0; j < ndim; ++j) {
    omega[j] = angVel*angNormal[j];
    alpha[j] = angAcc*angNormal[j];
    rotation[j] = 0.0;
  }

  for (int i = 0; i < boundarySize; i++) {
    int index = boundaryID[i]*ndim;
    for (int j = 0; j < ndim; ++j) {
      position[j] = coordinates[index+j]+displacements[index+j]-cm[j];
    }
    for (int j = 0; j < ndim; ++j) {
      for (int k = 0; k < ndim; ++k) {
        rotation[j] += rotMat[j][k]*position[k];
      }
      rotation[j] -= position[j];
    }
    crossProduct(omega, position, velRotation);
    crossProduct(omega, linVel, accCorioli);
    crossProduct(omega, velRotation, centrifugal);
    crossProduct(alpha, position, accRotation);
    for (int j = 0; j < ndim; ++j) {
      displacements[index+j] += (linDispl[j]-linDisplOld[j])+rotation[j];
      velocities[index+j] = linVel[j]+velRotation[j];
      // For energy computations
      accelerations[index+j] = linAcc[j]+2.0*accCorioli[j]+accRotation[j]+centrifugal[j];
    }
  }
  if (world_rank == 0) {
    FILE *datFile;
    datFile = fopen("motion.dat", "a");
    fprintf(datFile, "%11.3e  %11.3e  %11.3e  %11.3e  %11.3e  %11.3e  %11.3e  %11.3e  %11.3e  %11.3e  %11.3e  %11.3e  %11.3e\n", Time, linAcc[0], linAcc[1], linAcc[2], linVel[0], linVel[1], linVel[2], linDispl[0], linDispl[1], linDispl[2], angAcc, angVel, angDispl*radToDeg);
    fclose(datFile);
  }
  for (int j = 0; j < ndim; ++j) {
    linDisplOld[j] = linDispl[j];
  }
  return;
}

void InitCustomPlot() {
  double xPlot = 0.0;
  double yPlot = 0.0;
  double zPlot = 0.0;
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
  // TODO : If multiple points have same point to plot use the lowest rank
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

void InitBoundaryCondition(double *aMax, double angMax) {
  double tol = 1e-5;
  // Find count of nodes with specified partID
  int rigidNodeCount = 0;
  for (int i = 0; i < nelements; ++i) {
    if (pid[i] == rigidPartID) {
      rigidNodeCount = rigidNodeCount + (eptr[i+1]-eptr[i]);
    }
  }
  // Allocate node storage
  int *rigidNodeID = (int*)malloc(rigidNodeCount*sizeof(int));
  if (rigidNodeID == NULL) {
    printf("ERROR(%d) : Unable to alocate rigidNodeID\n", world_rank);
    exit(0);
  }
  // Store all nodes to be made rigid
  int nodePtr = 0;
  for (int i = 0; i < nelements; ++i) {
    if (pid[i] == rigidPartID) {
      for (int j = eptr[i]; j < eptr[i+1]; ++j) {
        rigidNodeID[nodePtr] = connectivity[j];
        nodePtr = nodePtr + 1;
      }
    }
  }
  assert(nodePtr == rigidNodeCount);
  // Sort and make unique
  qsort(rigidNodeID, rigidNodeCount, sizeof(int), compare);
  boundarySize = unique(rigidNodeID, rigidNodeCount);
  printf("INFO(%d): %d nodes given rigid motion\n", world_rank, boundarySize);
  boundaryID = (int*)malloc(boundarySize*sizeof(int));
  if (boundaryID == NULL) {
    printf("ERROR(%d) : Unable to alocate boundaryID\n", world_rank);
    exit(0);
  }
  for (int i = 0; i < boundarySize; ++i) {
    int node = rigidNodeID[i];
    boundaryID[i] = node;
    boundary[node] = 1;
  }
  free(rigidNodeID);
  // Compute the constants required for acceleration computations
  for (int i = 0; i < ndim; ++i) {
    aLin[i] = aMax[i]/peakTime;
    bLin[i] = aMax[i]/(tMax-peakTime);
  }
  aAng = angMax/peakTime;
  bAng = angMax/(tMax-peakTime);
  for (int i = 0; i < 3; ++i) {
    linDisplOld[i] = 0.0;
  }

  FILE *datFile;
  if (world_rank == 0) {
    datFile = fopen("motion.dat", "w");
    fprintf(datFile, "# Motion Results\n");
    fprintf(datFile, "# Time     AccX    AccY    AccZ    VelX    VelY    VelZ    LinDispX  LinDisp Y  LinDisp Z  AngAcc   AngVel   Angle \n");
    fclose(datFile);
  }

  return;
}
