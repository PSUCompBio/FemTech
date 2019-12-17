#include "FemTech.h"
#include "blas.h"
#include "gitbranch.h"
#include "jsonfuncs.h"

#include <assert.h>

/*Declare Functions*/
void CustomPlot();
void InitCustomPlot();
void InitBoundaryCondition(double *aMax, double angMax);
void updateBoundaryNeighbour(void);
void ApplyAccBoundaryConditions();
void WriteMaxStrainFile(double maxStrain, double maxX, double maxY, \
    double maxZ, double maxT, double minStrain, double minX, double minY, \
    double minZ, double minT);

/* Global Variables/Parameters */
double Time, dt;
int nSteps;
bool ImplicitStatic = false;
bool ImplicitDynamic = false;
bool ExplicitDynamic = true;
double ExplicitTimeStepReduction = 0.8;
double FailureTimeStep = 1e-11;
static const double radToDeg = 180.0 / (atan(1.0) * 4.0);
int nPlotSteps = 50;

/* Global variables used only in this file */
int nodeIDtoPlot;
bool rankForCustomPlot;
/* Global variables for bc */
int *boundaryID = NULL;
int boundarySize;
double aLin[3], bLin[3];
double aAng, bAng;
double peakTime, tMax;
double thetaOld = 0.0;
double linDisplOld[3];
double angNormal[3];
const int rigidPartID = 2; // part ID of elements to be made rigid

int main(int argc, char **argv) {
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);
  // Get the number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  // Get the rank of the process
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  static const double gC = 9.81;

  Json::Value simulationJson = getConfig(argv[1]);
  std::string meshFile = simulationJson["mesh"].asString();
  double accMax[3];
  accMax[0] = gC*simulationJson["linear-acceleration"][0].asDouble();
  accMax[1] = gC*simulationJson["linear-acceleration"][1].asDouble();
  accMax[2] = gC*simulationJson["linear-acceleration"][2].asDouble();

  double angAccMax = simulationJson["angular-acceleration"].asDouble();
  angNormal[0] = 0.0;
  angNormal[1] = 0.0;
  angNormal[2] = 1.0;
  peakTime = 0.020;
  tMax = 0.040;
  if (world_rank == 0) {
    printf("INFO : Git commit : %s of branch %s\n", GIT_COMMIT_HASH,
           GIT_BRANCH);
    printf("INFO : Linear Acceleration : (%7.3e, %7.3e, %7.3e)\n", accMax[0],
           accMax[1], accMax[2]);
    printf("INFO : Angular Acceleration : %7.3e\n", angAccMax);
    printf("INFO : Reading Mesh File : %s\n", meshFile.c_str());
  }

  ReadInputFile(meshFile.c_str());
  ReadMaterials();

  PartitionMesh();

  AllocateArrays();
  InitCustomPlot();
  InitBoundaryCondition(accMax, angAccMax);

  /* Write inital, undeformed configuration*/
  Time = 0.0;
  int plot_counter = 0;
  WriteVTU(meshFile.c_str(), plot_counter);
  stepTime[plot_counter] = Time;
  CustomPlot();

  // Dynamic Explcit solution using....
  dt = 0.0;

  int time_step_counter = 0;
  const int nDOF = nnodes * ndim;
  /** Central Difference Method - Beta and Gamma */
  // double beta = 0;
  // double gamma = 0.5;

  ShapeFunctions();
  /*  Step-1: Calculate the mass matrix similar to that of belytschko. */
  AssembleLumpedMass();

  // Used if initial velocity and acceleration BC is to be set.

  ApplyAccBoundaryConditions();

  /* Obtain dt, according to Belytschko dt is calculated at end of getForce */
  dt = ExplicitTimeStepReduction * StableTimeStep();
  /* Step-2: getforce step from Belytschko */
  GetForce(); // Calculating the force term.

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

  /* Variables to compute maximim and minimum strain */
  double maxStrain = 0.0, minStrain = 0.0;
  int maxElem = 0, minElem = 0;
  double maxT = 0.0, minT = 0.0;
  struct {
    double value;
    int   rank;
  } parStructMax, parStructMin; 

  /* Step-4: Time loop starts....*/
  while (Time < tMax) {
    t_n = Time;
    double t_np1 = Time + dt;
    Time = t_np1; /*Update the time by adding full time step */
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
    int writeFlag = time_step_counter%nsteps_plot;
    CheckEnergy(Time, writeFlag);

    if (writeFlag == 0) {
      if (world_rank == 0) {
        printf("Time : %15.6e, dt=%15.6e, tmax : %15.6e\n", Time, dt, tMax);
      }
      plot_counter = plot_counter + 1;
      double currentStrainMaxElem, currentStrainMinElem;
      double currentStrainMax = 0.0, currentStrainMin = 0.0;
      for (int i = 0; i < nelements; i++) {
        CalculateMaximumPrincipalStrain(i, &currentStrainMaxElem, &currentStrainMinElem);
        if (currentStrainMax < currentStrainMaxElem) {
          currentStrainMax = currentStrainMaxElem;
          maxElem = i;
        }
        if (currentStrainMin > currentStrainMinElem) {
          currentStrainMin = currentStrainMinElem;
          minElem = i;
        }
      } // calculating max and minimum strain over local elements
      // Updating max and min time
      if (currentStrainMax > maxStrain) {
        maxT = Time;
        maxStrain = currentStrainMax;
      }
      if (currentStrainMin < minStrain) {
        minT = Time;
        minStrain = currentStrainMin;
      }
      // Find the gloabl min and max strain
      parStructMax.value = maxStrain;
      parStructMax.rank = world_rank;
      MPI_Allreduce(MPI_IN_PLACE, &parStructMax, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
      maxStrain = parStructMax.value;

      parStructMin.value = minStrain;
      parStructMin.rank = world_rank;
      MPI_Allreduce(MPI_IN_PLACE, &parStructMin, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
      minStrain = parStructMin.value;

      printf("------Plot %d: WriteVTU by rank : %d\n", plot_counter,
             world_rank);
      WriteVTU(meshFile.c_str(), plot_counter);
      if (plot_counter < MAXPLOTSTEPS) {
        stepTime[plot_counter] = Time;
        WritePVD(meshFile.c_str(), plot_counter);
      }
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
  // Copy max and min strain 
  if (parStructMin.rank == world_rank) {
    double minLocationAndTime[4];
    for (int i = 0; i < 4; ++i) {
      minLocationAndTime[i] = 0.0;
    }
    // Compute element coordinates 
    int nP = eptr[minElem+1]-eptr[minElem];
    for (int i = eptr[minElem]; i < eptr[minElem+1]; ++i) {
      minLocationAndTime[0] += coordinates[connectivity[i]*ndim];
      minLocationAndTime[1] += coordinates[connectivity[i]*ndim+1];
      minLocationAndTime[2] += coordinates[connectivity[i]*ndim+2];
    }
    for (int i = 0; i < ndim; ++i) {
      minLocationAndTime[i] = minLocationAndTime[i]/((double)nP);
    }
    minLocationAndTime[3] = minT;
    MPI_Send(minLocationAndTime, 4, MPI_DOUBLE, 0, 7297, MPI_COMM_WORLD);
  }
  if (parStructMax.rank == world_rank) {
    double maxLocationAndTime[4];
    for (int i = 0; i < 4; ++i) {
      maxLocationAndTime[i] = 0.0;
    }
    // Compute element coordinates 
    int nP = eptr[maxElem+1]-eptr[maxElem];
    for (int i = eptr[maxElem]; i < eptr[maxElem+1]; ++i) {
      maxLocationAndTime[0] += coordinates[connectivity[i]*ndim];
      maxLocationAndTime[1] += coordinates[connectivity[i]*ndim+1];
      maxLocationAndTime[2] += coordinates[connectivity[i]*ndim+2];
    }
    for (int i = 0; i < ndim; ++i) {
      maxLocationAndTime[i] = maxLocationAndTime[i]/((double)nP);
    }
    maxLocationAndTime[3] = maxT;
    MPI_Send(maxLocationAndTime, 4, MPI_DOUBLE, 0, 7298, MPI_COMM_WORLD);
  }
  if (world_rank == 0) {
    double minRecv[4];
    double maxRecv[4];
    MPI_Recv(minRecv, 4, MPI_DOUBLE, parStructMin.rank, 7297, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(maxRecv, 4, MPI_DOUBLE, parStructMax.rank, 7298, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // if (plot_counter < MAXPLOTSTEPS) {
    //   stepTime[plot_counter] = Time;
    //   WritePVD(meshFile.c_str(), plot_counter, Time);
    // }
    WriteMaxStrainFile(maxStrain, maxRecv[0], maxRecv[1], maxRecv[2], \
        maxRecv[3], minStrain, minRecv[0], minRecv[1], minRecv[2], minRecv[3]);
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
      linAcc[i] = aLin[i] * Time;
      linVel[i] = 0.5 * aLin[i] * Time * Time;
      linDispl[i] = aLin[i] * Time * Time * Time / 6.0;
    }

    angAcc = aAng * Time;
    angVel = 0.5 * aAng * Time * Time;
    angDispl = aAng * Time * Time * Time / 6.0;
  } else {
    if (Time < tMax) {
      for (int i = 0; i < ndim; ++i) {
        linAcc[i] = (aLin[i] + bLin[i]) * peakTime - bLin[i] * Time;
        linVel[i] = (aLin[i] + bLin[i]) *
                        (peakTime * Time - 0.5 * peakTime * peakTime) -
                    0.5 * bLin[i] * Time * Time;
        linDispl[i] =
            0.5 * (aLin[i] + bLin[i]) * peakTime *
                (peakTime * peakTime / 3.0 - peakTime * Time + Time * Time) -
            bLin[i] * Time * Time * Time / 6.0;
      }

      angAcc = (aAng + bAng) * peakTime - bAng * Time;
      angVel = (aAng + bAng) * (peakTime * Time - 0.5 * peakTime * peakTime) -
               0.5 * bAng * Time * Time;
      angDispl =
          0.5 * (aAng + bAng) * peakTime *
              (peakTime * peakTime / 3.0 - peakTime * Time + Time * Time) -
          bAng * Time * Time * Time / 6.0;
    } else {
      for (int i = 0; i < ndim; ++i) {
        linAcc[i] = 0.0;
        linVel[i] = (aLin[i] + bLin[i]) *
                        (peakTime * tMax - 0.5 * peakTime * peakTime) -
                    0.5 * bLin[i] * tMax * tMax;
        linDispl[i] =
            0.5 * (aLin[i] + bLin[i]) * peakTime *
                (peakTime * peakTime / 3.0 - peakTime * tMax + tMax * tMax) -
            bLin[i] * tMax * tMax * tMax / 6.0 + linVel[i] * (Time - tMax);
      }
      angAcc = 0.0;
      angVel = (aAng + bAng) * (peakTime * tMax - 0.5 * peakTime * peakTime) -
               0.5 * bAng * tMax * tMax;
      angDispl =
          0.5 * (aAng + bAng) * peakTime *
              (peakTime * peakTime / 3.0 - peakTime * tMax + tMax * tMax) -
          bAng * tMax * tMax * tMax / 6.0 + angVel * (Time - tMax);
    }
  }
  double dTheta = angDispl - thetaOld;
  thetaOld = angDispl;
  GetBodyCenterofMass(cm);
  get3dRotationMatrix(angNormal, dTheta, rotMat);
  for (int j = 0; j < ndim; ++j) {
    omega[j] = angVel * angNormal[j];
    alpha[j] = angAcc * angNormal[j];
  }

  for (int i = 0; i < boundarySize; i++) {
    int index = boundaryID[i] * ndim;
    for (int j = 0; j < ndim; ++j) {
      position[j] = coordinates[index + j] + displacements[index + j] - cm[j];
    }
    for (int j = 0; j < ndim; ++j) {
      rotation[j] = 0.0;
    }
    for (int j = 0; j < ndim; ++j) {
      for (int k = 0; k < ndim; ++k) {
        rotation[j] += rotMat[j][k] * position[k];
      }
      rotation[j] -= position[j];
    }
    crossProduct(omega, position, velRotation);
    crossProduct(omega, linVel, accCorioli);
    crossProduct(omega, velRotation, centrifugal);
    crossProduct(alpha, position, accRotation);
    for (int j = 0; j < ndim; ++j) {
      displacements[index + j] += (linDispl[j] - linDisplOld[j]) + rotation[j];
      velocities[index + j] = linVel[j] + velRotation[j];
      // For energy computations
      accelerations[index + j] =
          linAcc[j] + 2.0 * accCorioli[j] + accRotation[j] + centrifugal[j];
    }
  }
  if (world_rank == 0) {
    FILE *datFile;
    datFile = fopen("motion.dat", "a");
    fprintf(datFile,
            "%11.3e  %11.3e  %11.3e  %11.3e  %11.3e  %11.3e  %11.3e  %11.3e  "
            "%11.3e  %11.3e  %11.3e  %11.3e  %11.3e\n",
            Time, linAcc[0], linAcc[1], linAcc[2], linVel[0], linVel[1],
            linVel[2], linDispl[0], linDispl[1], linDispl[2], angAcc, angVel,
            angDispl * radToDeg);
    fclose(datFile);
  }
  for (int j = 0; j < ndim; ++j) {
    linDisplOld[j] = linDispl[j];
  }
  return;
}

void InitCustomPlot() {
  double xPlot = 2.087348700e-02;
  double yPlot = 2.087348700e-02;
  double zPlot = 2.175215120e-02;
  double tol = 1e-5;

  int idToPlot;
  FILE *datFile;
  rankForCustomPlot = false;
  int index;
  const int x = 0;
  const int y = 1;
  const int z = 2;
  int state = 0;
  int cumState = 0;

  for (int i = 0; i < nnodes; ++i) {
    if (fabs(coordinates[ndim * i + x] - xPlot) < tol &&
        fabs(coordinates[ndim * i + y] - yPlot) < tol &&
        fabs(coordinates[ndim * i + z] - zPlot) < tol) {
      nodeIDtoPlot = i;
      rankForCustomPlot = true;
      state = 1;
      break;
    }
  }
  // If multiple procs have same node use the lowest rank
  MPI_Scan(&state, &cumState, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (rankForCustomPlot) {
    if (cumState > 1) {
      rankForCustomPlot = false;
      return;
    }
    printf("INFO(%d) : nodeID for plot : %d (%15.9e, %15.9e, %15.9e)\n",
          world_rank, nodeIDtoPlot, coordinates[ndim * nodeIDtoPlot + x],
          coordinates[ndim * nodeIDtoPlot + y],
          coordinates[ndim * nodeIDtoPlot + z]);
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
    const int index = nodeIDtoPlot * ndim;

    FILE *datFile;
    datFile = fopen("plot.dat", "a");
    fprintf(datFile, "%15.9e %15.9e  %15.9e  %15.9e\n", Time,
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
      rigidNodeCount = rigidNodeCount + (eptr[i + 1] - eptr[i]);
    }
  }
  // Allocate node storage
  int *rigidNodeID = (int *)malloc(rigidNodeCount * sizeof(int));
  if (rigidNodeID == NULL) {
    printf("ERROR(%d) : Unable to alocate rigidNodeID\n", world_rank);
    exit(0);
  }
  // Store all nodes to be made rigid
  int nodePtr = 0;
  for (int i = 0; i < nelements; ++i) {
    if (pid[i] == rigidPartID) {
      for (int j = eptr[i]; j < eptr[i + 1]; ++j) {
        rigidNodeID[nodePtr] = connectivity[j];
        nodePtr = nodePtr + 1;
      }
    }
  }
  assert(nodePtr == rigidNodeCount);
  // Sort and make unique
  qsort(rigidNodeID, rigidNodeCount, sizeof(int), compare);
  boundarySize = unique(rigidNodeID, rigidNodeCount);
  for (int i = 0; i < boundarySize; ++i) {
    int index = rigidNodeID[i] * ndim;
    boundary[index] = 1;
    boundary[index + 1] = 1;
    boundary[index + 2] = 1;
  }
  updateBoundaryNeighbour();
  free(rigidNodeID);
  boundarySize = 0;
  for (int i = 0; i < nnodes*ndim; i += 3) {
    if (boundary[i]) {
      boundarySize = boundarySize + 1;
    }
  }

  printf("INFO(%d): %d nodes given rigid motion\n", world_rank, boundarySize);
  boundaryID = (int *)malloc(boundarySize * sizeof(int));
  if (boundaryID == NULL) {
    printf("ERROR(%d) : Unable to alocate boundaryID\n", world_rank);
    exit(0);
  }

  int idIndex = 0;
  for (int i = 0; i < nnodes; ++i) {
    int index = i * ndim;
    if (boundary[index]) {
      boundaryID[idIndex] = i;
      idIndex = idIndex + 1;
    }
  }
  assert(idIndex == boundarySize);
  // Compute the constants required for acceleration computations
  for (int i = 0; i < ndim; ++i) {
    aLin[i] = aMax[i] / peakTime;
    bLin[i] = aMax[i] / (tMax - peakTime);
  }
  aAng = angMax / peakTime;
  bAng = angMax / (tMax - peakTime);
  for (int i = 0; i < 3; ++i) {
    linDisplOld[i] = 0.0;
  }

  FILE *datFile;
  if (world_rank == 0) {
    datFile = fopen("motion.dat", "w");
    fprintf(datFile, "# Motion Results\n");
    fprintf(datFile,
            "# Time     AccX    AccY    AccZ    VelX    VelY    VelZ    "
            "LinDispX  LinDisp Y  LinDisp Z  AngAcc   AngVel   Angle \n");
    fclose(datFile);
  }
  return;
}

void WriteMaxStrainFile(double maxStrain, double maxX, double maxY, \
    double maxZ, double maxT, double minStrain, double minX, double minY, \
    double minZ, double minT) {
  // Write the maximum strain to file
  FILE *maxStrainFile;
  maxStrainFile = fopen("maxstrain.dat", "w");
  fprintf(maxStrainFile, "# Maximum and Minimum Principal Strain Results\n");
  fprintf(maxStrainFile, "# MaxStrain  MaxX  MaxY MaxZ  MaxT  MinStrain  MinX  MinY MinZ  MinT\n");
  fprintf(maxStrainFile,
          "%11.6e  %11.6e  %11.6e  %11.6e  %11.6e  %11.6e  %11.6e  %11.6e  "
          "%11.6e  %11.6e\n", maxStrain, maxX, maxY, maxZ, maxT, minStrain, minX, minY, minZ, minT);
  fclose(maxStrainFile);
}

void updateBoundaryNeighbour(void) {
  int* sendNodeBoundary;
  int* recvNodeBoundary;
  int totalNodeToSend = sendNeighbourCountCum[sendProcessCount];
  recvNodeBoundary = (int*)malloc(ndim*totalNodeToSend * sizeof(int));
  sendNodeBoundary = (int*)malloc(ndim*totalNodeToSend * sizeof(int));
  // Update array to send 
  for (int i = 0; i < totalNodeToSend; ++i) {
    int nodeIndex = ndim*sendNodeIndex[i];
    sendNodeBoundary[ndim*i] = boundary[nodeIndex];
    sendNodeBoundary[ndim*i+1] = boundary[nodeIndex+1];
    sendNodeBoundary[ndim*i+2] = boundary[nodeIndex+2];
  }
  // Create send requests
  int boundaryTag = 2798;
  MPI_Request *requestListSend =
        (MPI_Request *)malloc(sizeof(MPI_Request) * sendProcessCount);
  for (int i = 0; i < sendProcessCount; ++i) {
    int process = sendProcessID[i];
    int location = sendNeighbourCountCum[i]*ndim;
    int size = ndim*sendNeighbourCount[i];
    MPI_Isend(&(sendNodeBoundary[location]), size, MPI_INT, process, \
        boundaryTag, MPI_COMM_WORLD, &(requestListSend[i]));
  }
  // Create recv requests
  MPI_Request *requestListRecv =
        (MPI_Request *)malloc(sizeof(MPI_Request) * sendProcessCount);
  for (int i = 0; i < sendProcessCount; ++i) {
    int process = sendProcessID[i];
    int location = sendNeighbourCountCum[i]*ndim;
    int size = ndim*sendNeighbourCount[i];
    MPI_Irecv(&(recvNodeBoundary[location]), size, MPI_INT, process, \
        boundaryTag, MPI_COMM_WORLD, &(requestListRecv[i]));
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
    if (recvNodeBoundary[ndim*i]) {
      boundary[nodeIndex] = 1;
    }
    if (recvNodeBoundary[ndim*i+1]) {
      boundary[nodeIndex+1] = 1;
    }
    if (recvNodeBoundary[ndim*i+2]) {
      boundary[nodeIndex+2] = 1;
    }
  }
  free(requestListSend);
  free(requestListRecv);
  free(recvNodeBoundary);
  free(sendNodeBoundary);
}
