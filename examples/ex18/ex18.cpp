#include "FemTech.h"
#include "blas.h"
#include "gitbranch.h"
#include "jsonfuncs.h"
#include "utilities.h"

#include <assert.h>

/*Declare Functions*/
void CustomPlot();
void InitCustomPlot();
void InitBoundaryCondition(const Json::Value& jsonInput);
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
double peakTime, tMax;
double thetaOld = 0.0;
double linDisplOld[3];
const int rigidPartID = 0; // part ID of elements to be made rigid

/* Variables used to store acceleration values */
int linAccXSize, linAccYSize, linAccZSize;
int angAccXSize, angAccYSize, angAccZSize;
double *linAccXt, *linAccYt, *linAccZt;
double *linAccXv, *linAccYv, *linAccZv;
double *angAccXt, *angAccYt, *angAccZt;
double *angAccXv, *angAccYv, *angAccZv;

int main(int argc, char **argv) {
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);
  // Get the number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  // Get the rank of the process
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // Read the input file
  Json::Value simulationJson = getConfig(argv[1]);
  std::string meshFile = simulationJson["mesh"].asString();

  peakTime = simulationJson["time-peak-acceleration"].asDouble();
  tMax = simulationJson["maximum-time"].asDouble();
  if (world_rank == 0) {
    printf("INFO : Git commit : %s of branch %s\n", GIT_COMMIT_HASH,
           GIT_BRANCH);
    printf("INFO : Reading Mesh File : %s\n", meshFile.c_str());
  }

  ReadInputFile(meshFile.c_str());
  // Read material properties before mesh partition to estimate material type kernel compute intensity
  ReadMaterials();

  PartitionMesh();

  AllocateArrays();
  InitCustomPlot();
  InitBoundaryCondition(simulationJson);

  /* Write inital, undeformed configuration*/
  Time = 0.0;
  int plot_counter = 0;
  WriteVTU(meshFile.c_str(), plot_counter);
  stepTime[plot_counter] = Time;
  CustomPlot();

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
  // In GetForce for viscoelastic material dt is required, hence we compute dt
  // prior to getforce to avoid special treatment of getforce at Time = 0
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
  // Free local boundary condition related arrays
  free1DArray(boundaryID);
  free1DArray(linAccXt);
  free1DArray(linAccYt);
  free1DArray(linAccZt);
  free1DArray(linAccXv);
  free1DArray(linAccYv);
  free1DArray(linAccZv);
  free1DArray(angAccXt);
  free1DArray(angAccYt);
  free1DArray(angAccZt);
  free1DArray(angAccXv);
  free1DArray(angAccYv);
  free1DArray(angAccZv);
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

  double dTheta = angDispl - thetaOld;
  thetaOld = angDispl;
  GetBodyCenterofMass(cm);

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
  double xPlot = 0.0;
  double yPlot = -0.366716;
  double zPlot = -0.018784;
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

void InitBoundaryCondition(const Json::Value& jsonInput) {
  // Read input JSON for acceleration values
  // Read linear acceleration and angular acceleration time traces
  linAccXSize = jsonInput["linear-acceleration"]["xt"].size();
  int tempSize = jsonInput["linear-acceleration"]["xv"].size();
  assert(tempSize == linAccXSize);
  linAccXt = (double*)malloc(sizeof(double)*linAccXSize);
  linAccXv = (double*)malloc(sizeof(double)*linAccXSize);
  linAccYSize = jsonInput["linear-acceleration"]["yt"].size();
  tempSize = jsonInput["linear-acceleration"]["yv"].size();
  assert(tempSize == linAccYSize);
  linAccYt = (double*)malloc(sizeof(double)*linAccYSize);
  linAccYv = (double*)malloc(sizeof(double)*linAccYSize);
  linAccZSize = jsonInput["linear-acceleration"]["zt"].size();
  tempSize = jsonInput["linear-acceleration"]["zv"].size();
  assert(tempSize == linAccYSize);
  linAccZt = (double*)malloc(sizeof(double)*linAccZSize);
  linAccZv = (double*)malloc(sizeof(double)*linAccZSize);
  jsonToArray(linAccXt, jsonInput["linear-acceleration"]["xt"]);
  jsonToArray(linAccXv, jsonInput["linear-acceleration"]["xv"]);
  jsonToArray(linAccYt, jsonInput["linear-acceleration"]["yt"]);
  jsonToArray(linAccYv, jsonInput["linear-acceleration"]["yv"]);
  jsonToArray(linAccZt, jsonInput["linear-acceleration"]["zt"]);
  jsonToArray(linAccZv, jsonInput["linear-acceleration"]["zv"]);

  angAccXSize = jsonInput["angular-acceleration"]["xt"].size();
  tempSize = jsonInput["angular-acceleration"]["xv"].size();
  assert(tempSize == angAccXSize);
  angAccXt = (double*)malloc(sizeof(double)*angAccXSize);
  angAccXv = (double*)malloc(sizeof(double)*angAccXSize);
  angAccYSize = jsonInput["angular-acceleration"]["yt"].size();
  tempSize = jsonInput["angular-acceleration"]["yv"].size();
  assert(tempSize == angAccYSize);
  angAccYt = (double*)malloc(sizeof(double)*angAccYSize);
  angAccYv = (double*)malloc(sizeof(double)*angAccYSize);
  angAccZSize = jsonInput["angular-acceleration"]["zt"].size();
  tempSize = jsonInput["angular-acceleration"]["zv"].size();
  assert(tempSize == angAccZSize);
  angAccZt = (double*)malloc(sizeof(double)*angAccZSize);
  angAccZv = (double*)malloc(sizeof(double)*angAccZSize);
  jsonToArray(angAccXt, jsonInput["angular-acceleration"]["xt"]);
  jsonToArray(angAccXv, jsonInput["angular-acceleration"]["xv"]);
  jsonToArray(angAccYt, jsonInput["angular-acceleration"]["yt"]);
  jsonToArray(angAccYv, jsonInput["angular-acceleration"]["yv"]);
  jsonToArray(angAccZt, jsonInput["angular-acceleration"]["zt"]);
  jsonToArray(angAccZv, jsonInput["angular-acceleration"]["zv"]);
  printf("Size of Xt = %d, lastValue = %f, %f\n", linAccXSize, linAccXt[linAccXSize-1], linAccXv[linAccXSize-1]);
  printf("Size of Yt = %d\n", linAccYSize);
  printf("Size of Zt = %d\n", linAccZSize);
  printf("Size of Xt = %d\n", angAccXSize);
  printf("Size of Yt = %d\n", angAccYSize);
  printf("Size of Zt = %d\n", angAccZSize);
  exit(0);

  static const double gC = 9.81;
  // Convert linear accelerations from g force to m/s^2
  // Convert time from milli-seconds to seconds
  for (int i = 0; i < linAccXSize; ++i) {
    linAccXv[i] = gC*linAccXv[i];
    linAccXt[i] = 0.001*linAccXt[i];
  }
  for (int i = 0; i < linAccYSize; ++i) {
    linAccYv[i] = gC*linAccYv[i];
    linAccYt[i] = 0.001*linAccYt[i];
  }
  for (int i = 0; i < linAccZSize; ++i) {
    linAccZv[i] = gC*linAccZv[i];
    linAccZt[i] = 0.001*linAccZt[i];
  }
  for (int i = 0; i < angAccXSize; ++i) {
    angAccXt[i] = 0.001*angAccXt[i];
  }
  for (int i = 0; i < angAccYSize; ++i) {
    angAccYt[i] = 0.001*angAccYt[i];
  }
  for (int i = 0; i < angAccZSize; ++i) {
    angAccZt[i] = 0.001*angAccZt[i];
  }

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
