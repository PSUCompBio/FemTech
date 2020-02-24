#include "FemTech.h"
#include "jsonfuncs.h"
#include "utilities.h"

#include "json/writer.h"

#include <map>
#include <string>
#include <fstream> // For JSON output
#include <memory> // For JSON output unique pointer
#include <math.h>
#include <assert.h>

#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

/*Declare Functions*/
void CustomPlot();
void InitCustomPlot();
void InitBoundaryCondition(const Json::Value& jsonInput);
void updateBoundaryNeighbour(void);
void ApplyAccBoundaryConditions();
int getImpactID(std::string location);
void WriteOutputFile(std::string);
void CalculateInjuryCriterions(void);

/* Global Variables/Parameters */
double Time, dt;
int nSteps;
bool ImplicitStatic = false;
bool ImplicitDynamic = false;
bool ExplicitDynamic = true;
double ExplicitTimeStepReduction = 0.8;
double FailureTimeStep = 5e-8; // Set for max runtime of around 5 hrs on aws
int nPlotSteps = 50;

/* Global variables used only in this file */
int nodeIDtoPlot;
bool rankForCustomPlot;
/* Global variables for bc */
int *boundaryID = NULL;
int boundarySize;
double peakTime, tMax;

/* Variables to compute maximim and minimum strain */
double maxStrain = 0.0, minStrain = 0.0, maxShear = 0.0;
int maxElem = 0, minElem = 0, shearElem = 0;
double maxT = 0.0, minT = 0.0, maxShearT = 0.0;
struct {
  double value;
  int   rank;
} parStructMax, parStructMin, parStructSMax; 

/* Variables used to store acceleration values */
int linAccXSize, linAccYSize, linAccZSize;
int angAccXSize, angAccYSize, angAccZSize;
double *linAccXt, *linAccYt, *linAccZt;
double *linAccXv, *linAccYv, *linAccZv;
double *angAccXt, *angAccYt, *angAccZt;
double *angAccXv, *angAccYv, *angAccZv;

/* Varibles for Local integrator */
const int nIntVar = 12;
typedef boost::array<double, nIntVar> state_type;
runge_kutta_dopri5<state_type> rk;
void computeDerivatives(const state_type &y, state_type &ydot, const double t);
state_type yInt, ydotInt;
double cm[3];

int main(int argc, char **argv) {
  // Initialize FemTech including logfile and MPI
  std::string uid = InitFemTech(argc, argv);
  // Read the input file
  Json::Value simulationJson = getConfig(argv[1]);

  std::string meshFile = simulationJson["mesh"].asString();
  tMax = simulationJson["maximum-time"].asDouble();
  FILE_LOG_MASTER(INFO, "Reading Mesh File : %s", meshFile.c_str());
  // Read Input Mesh file and equally partition elements among processes
  ReadInputFile(meshFile.c_str());
  size_t lastindex = meshFile.find_last_of(".");
  std::string outputFileName = meshFile.substr(0, lastindex) + "_" + uid;
  // Read material properties before mesh partition to estimate 
  // material type kernel compute intensity
  ReadMaterials();

  PartitionMesh();

  AllocateArrays();
  // InitCustomPlot();

  // Initial settings for BC evaluations
  // Used if initial velocity and acceleration BC is to be set.
  InitBoundaryCondition(simulationJson);

  /* Write inital, undeformed configuration*/
  Time = 0.0;
  int plot_counter = 0;
  WriteVTU(outputFileName.c_str(), plot_counter);
  stepTime[plot_counter] = Time;
  // CustomPlot();

  int time_step_counter = 0;
  /** Central Difference Method - Beta and Gamma */
  // double beta = 0;
  // double gamma = 0.5;

  ShapeFunctions();
  /*  Step-1: Calculate the mass matrix similar to that of belytschko. */
  AssembleLumpedMass();

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

  FILE_LOG_MASTER(INFO, "initial dt = %3.3e, nSteps = %d, nsteps_plot = %d", dt, nSteps,
           nsteps_plot);

  time_step_counter = time_step_counter + 1;
  double t_n = 0.0;

  FILE_LOG_MASTER(INFO, "------------------------------- Loop ----------------------------");
  FILE_LOG_MASTER(INFO, "Time : %15.6e, tmax : %15.6e", Time, tMax);

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
      FILE_LOG_MASTER(INFO, "Time : %15.6e, dt=%15.6e, tmax : %15.6e", Time, dt, tMax);
      plot_counter = plot_counter + 1;
      CalculateInjuryCriterions();

      FILE_LOG(INFO, "------ Plot %d: WriteVTU", plot_counter);
      WriteVTU(outputFileName.c_str(), plot_counter);
      if (plot_counter < MAXPLOTSTEPS) {
        stepTime[plot_counter] = Time;
        WritePVD(outputFileName.c_str(), plot_counter);
      }
      // CustomPlot();
      FILE_LOGMatrixRM(DEBUGLOG, displacements, nNodes, ndim, "Displacement Solution");
    }
    time_step_counter = time_step_counter + 1;
    dt = ExplicitTimeStepReduction * StableTimeStep();

    // Write out the last time step
    // CustomPlot();
    // Barrier not a must
    MPI_Barrier(MPI_COMM_WORLD);
  } // end explcit while loop
  FILE_LOGMatrixRM(DEBUGLOG, displacements, nNodes, ndim, "Final Displacement Solution");

  WriteOutputFile(uid);
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

  FinalizeFemTech();
  return 0;
}

void ApplyAccBoundaryConditions() {
  double r[4], R[4], Rinv[4], V[4], Vp[4]; //quaternions
  double omegaR[3], omega[3], omegaOmegaR[3], omegaVel[3], vel[3]; //vectors
  double alpha[3], alphaR[3], locV[3];

  rk.do_step(computeDerivatives, yInt, ydotInt, Time-dt, dt);
  r[0] = 0.0; r[1] = yInt[3]; r[2] = yInt[4]; r[3] = yInt[5];
  quaternionExp(r, R); // R = exp(r)
  quaternionInverse(R, Rinv);
  omega[0] = yInt[0]; omega[1] = yInt[1]; omega[2] = yInt[2];
  alpha[0] = ydotInt[0]; alpha[1] = ydotInt[1]; alpha[2] = ydotInt[2];
  vel[0] = yInt[6]; vel[1] = yInt[7]; vel[2] = yInt[8];

  for (int i = 0; i < boundarySize; i++) {
    int index = boundaryID[i] * ndim;
    for (int j = 0; j < ndim; ++j) {
      locV[j] = coordinates[index + j] - cm[j];
    }
    V[0] = 0.0; V[1] = locV[0]; V[2] = locV[1]; V[3] = locV[2];
    quaternionRotate(V, R, Rinv, Vp); // Vp = RVR^{-1}
    crossProduct(omega, &(Vp[1]), omegaR);
    crossProduct(omega, omegaR, omegaOmegaR);
    crossProduct(omega, vel, omegaVel);
    crossProduct(alpha, &(Vp[1]), alphaR);
    for (int j = 0; j < ndim; ++j) {
      // Displacement
      displacements[index+j] = Vp[j+1] - locV[j] + yInt[9+j];
      // Velocity
      velocities[index+j] = omegaR[j] + yInt[6+j];
      // Acceleration
      accelerations[index+j] = 2.0*omegaVel[j] + omegaOmegaR[j] + ydotInt[6+j] + alphaR[j];
    }
  }
  return;
}

void InitCustomPlot() {
  double xPlot = -0.009213;
  double yPlot = 0.046231;
  double zPlot = 0.007533;
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

  for (int i = 0; i < nNodes; ++i) {
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
    FILE_LOG_SINGLE(INFO, "nodeID for plot : %d (%15.9e, %15.9e, %15.9e)",
          nodeIDtoPlot, coordinates[ndim * nodeIDtoPlot + x],
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
  static const double gC = 9.81;
  GetBodyCenterofMass(cm);
  // Read input JSON for acceleration values
  if (jsonInput["linear-acceleration"].isObject()) { // Use time traces from file
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
  } else { // Read maximum values from input file
    double accMax[3];
    accMax[0] = gC*jsonInput["linear-acceleration"][0].asDouble();
    accMax[1] = gC*jsonInput["linear-acceleration"][1].asDouble();
    accMax[2] = gC*jsonInput["linear-acceleration"][2].asDouble();
    peakTime = jsonInput["time-peak-acceleration"].asDouble();
    linAccXSize = 3; linAccYSize = 3; linAccZSize = 3;
    linAccXt = (double*)calloc(linAccXSize, sizeof(double));
    linAccXv = (double*)calloc(linAccXSize, sizeof(double));
    linAccYt = (double*)calloc(linAccYSize, sizeof(double));
    linAccYv = (double*)calloc(linAccYSize, sizeof(double));
    linAccZt = (double*)calloc(linAccZSize, sizeof(double));
    linAccZv = (double*)calloc(linAccZSize, sizeof(double));
    linAccXt[1] = linAccYt[1] = linAccZt[1] = peakTime;
    linAccXt[2] = linAccYt[2] = linAccZt[2] = tMax;
    linAccXv[1] = accMax[0]; 
    linAccYv[1] = accMax[1]; 
    linAccZv[1] = accMax[2];
    if (jsonInput["angular-acceleration"].size() == 3) {
      accMax[0] = jsonInput["angular-acceleration"][0].asDouble();
      accMax[1] = jsonInput["angular-acceleration"][1].asDouble();
      accMax[2] = jsonInput["angular-acceleration"][2].asDouble();
    } else {
      double angNormal[3];
      int impactPointID = getImpactID(jsonInput["impact-point"].asString());
      // Compute the axis of rotation based on center of mass and impact point
      double impactNodeCoord[3];
      // Find if global node ID is present on the current process
      int nodeStatus = coordinateFromGlobalID(globalNodeID, impactPointID, nNodes, \
          impactNodeCoord);
      // Find lowest rank process with node ID
      int *nodeWithID = (int*)malloc(world_size*sizeof(int));
      MPI_Allgather(&nodeStatus, 1, MPI_INT, nodeWithID, 1, MPI_INT, MPI_COMM_WORLD);
      int nodeIDGlobal = -1;
      for (int i = 0; i < world_size; ++i) {
        if (nodeWithID[i]) {
          nodeIDGlobal = i;
          break;
        }
      }
      // Recieve node co-ordinates
      MPI_Bcast(impactNodeCoord, ndim, MPI_DOUBLE, nodeIDGlobal, MPI_COMM_WORLD);
      FILE_LOG_MASTER(INFO, "NodeID of impact : %d (%15.9e, %15.9e, %15.9e)", \
          impactPointID, impactNodeCoord[0], impactNodeCoord[1], \
          impactNodeCoord[2]);
      //Compute the axis of rotation
      double norm = 0.0;
      for (int i = 0; i < ndim; ++i) {
        double ds = impactNodeCoord[i]-cm[i];
        angNormal[i] = ds;
        norm += ds*ds;
      }
      norm = sqrt(norm);
      for (int i = 0; i < ndim; ++i) {
        angNormal[i] /= norm;
      }
      double angAccMax = jsonInput["angular-acceleration"].asDouble();
      accMax[0] = angAccMax*angNormal[0];
      accMax[1] = angAccMax*angNormal[1];
      accMax[2] = angAccMax*angNormal[2];
      free(nodeWithID);
    }
    angAccXSize = 3; angAccYSize = 3; angAccZSize = 3;
    angAccXt = (double*)calloc(angAccXSize, sizeof(double));
    angAccXv = (double*)calloc(angAccXSize, sizeof(double));
    angAccYt = (double*)calloc(angAccYSize, sizeof(double));
    angAccYv = (double*)calloc(angAccYSize, sizeof(double));
    angAccZt = (double*)calloc(angAccZSize, sizeof(double));
    angAccZv = (double*)calloc(angAccZSize, sizeof(double));
    angAccXt[1] = angAccYt[1] = angAccZt[1] = peakTime;
    angAccXt[2] = angAccYt[2] = angAccZt[2] = tMax;
    angAccXv[1] = accMax[0]; 
    angAccYv[1] = accMax[1]; 
    angAccZv[1] = accMax[2];
  }

  double tol = 1e-5;
  // Find count of nodes with specified partID
  int rigidNodeCount = 0;
  for (int i = 0; i < nelements; ++i) {
    if (materialID[pid[i]] == 0) {
      rigidNodeCount = rigidNodeCount + (eptr[i + 1] - eptr[i]);
    }
  }
  // Allocate node storage
  int *rigidNodeID = (int *)malloc(rigidNodeCount * sizeof(int));
  if (rigidNodeID == NULL) {
    FILE_LOG_SINGLE(ERROR, "Unable to alocate rigidNodeID");
    TerminateFemTech(12);
  }
  // Store all nodes to be made rigid
  int nodePtr = 0;
  for (int i = 0; i < nelements; ++i) {
    if (materialID[pid[i]] == 0) {
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
  for (int i = 0; i < nDOF; i += 3) {
    if (boundary[i]) {
      boundarySize = boundarySize + 1;
    }
  }

  FILE_LOG(INFO, "%d nodes given rigid motion", boundarySize);
  boundaryID = (int *)malloc(boundarySize * sizeof(int));
  if (boundaryID == NULL) {
    FILE_LOG_SINGLE(ERROR, "Unable to alocate boundaryID");
    TerminateFemTech(12);
  }

  int idIndex = 0;
  for (int i = 0; i < nNodes; ++i) {
    int index = i * ndim;
    if (boundary[index]) {
      boundaryID[idIndex] = i;
      idIndex = idIndex + 1;
    }
  }
  assert(idIndex == boundarySize);
  // Setting initial conditions
  for (int j = 0; j < nIntVar; ++j) {
    yInt[j] = 0.0;
    ydotInt[j] = 0.0;
  }
  // Set initial displacements, accelerations, velocities to zero
  for (int i = 0; i < boundarySize; i++) {
    int index = boundaryID[i] * ndim;
    for (int j = 0; j < ndim; ++j) {
      displacements[index+j] = 0.0;
      velocities[index+j] = 0.0;
      accelerations[index+j] = 0.0;
    }
  }
  return;
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

int getImpactID(std::string location) {
  const std::map<std::string, int> pointMap{ {"top-right", 5328}, \
    {"top-left", 1749}, {"front-low", 2575}, {"front-high", 1967}, \
    {"right-high", 5720}, {"bottom-front", 2575}, {"top-front", 1967}, \
    {"top-rear", 1873}};
    if (pointMap.find(location) == pointMap.end() ) {
      FILE_LOG_SINGLE(ERROR, "Impact location not found, check value of impact-point");
      TerminateFemTech(11);
    }
    return pointMap.at(location);
}

void computeDerivatives(const state_type &y, state_type &ydot, const double t) {
  // ydot[0-2] : Store angular acceleration, alpha
  // ydot[3-5] : Store derivative of rotation quaternion generator, rdot
  // ydot[6-8] : Linear acceleration of the center of mass, acc
  // ydot[9-11] : Linear velocity of the center of mass, vel
  // Storing alpha
  ydot[0] = interpolateLinear(angAccXSize, angAccXt, angAccXv, t);
  ydot[1] = interpolateLinear(angAccYSize, angAccYt, angAccYv, t);
  ydot[2] = interpolateLinear(angAccZSize, angAccZt, angAccZv, t);
  // Linear Acceleration
  ydot[6] = interpolateLinear(linAccXSize, linAccXt, linAccXv, t);
  ydot[7] = interpolateLinear(linAccYSize, linAccYt, linAccYv, t);
  ydot[8] = interpolateLinear(linAccZSize, linAccZt, linAccZv, t);
  // Linear Velocity
  ydot[9] = y[6];
  ydot[10] = y[7];
  ydot[11] = y[8];
  // rdot
  // Ref: The integration of angular velocity, Michael Boyle, 2017
  // https://arxiv.org/pdf/1604.08139.pdf
  double r[3];
  r[0] = y[3]; r[1] = y[4]; r[2] = y[5];
  double *rdot = &(ydot[3]);
  double rMagnitude = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
  if (rMagnitude < 1e-10) {
    // rdot = 0.5*omega
    rdot[0] = 0.5*y[0]; 
    rdot[1] = 0.5*y[1]; 
    rdot[2] = 0.5*y[2]; 
  } else {
    double rCotR = rMagnitude/tan(rMagnitude);
    double omega[3]; omega[0] = y[0]; omega[1] = y[1]; omega[2] = y[2];
    crossProduct(omega, r, rdot); // compute rdot = omega x r
    // Normalize r by its magnitude
    for (int i = 0; i < 3; ++i) {
      r[i] = r[i]/rMagnitude;
    }
    double rDotOmega = dotProduct3D(r, omega); // compure r.omega/||r||
    for (int i = 0; i < 3; ++i) {
      // Compute 0.5 omegaxr + omega*(0.5||r|| cot(||r||) + (1-||r||cot(||r||))
      // *(r.omega)*r/(2*||r||*||r||)
      rdot[i] = 0.5*(rdot[i] + rCotR*y[i] + (1.0-rCotR)*rDotOmega*r[i]);
    }
  }
}

void WriteOutputFile(std::string uid) {
  /* Calculate min and max strain location and send to master */
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
    MPI_Send(&global_eid[minElem], 1, MPI_INT, 0, 7297, MPI_COMM_WORLD);
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
    MPI_Send(&global_eid[maxElem], 1, MPI_INT, 0, 7298, MPI_COMM_WORLD);
  }
  if (parStructSMax.rank == world_rank) {
    double maxLocationAndTime[4];
    for (int i = 0; i < 4; ++i) {
      maxLocationAndTime[i] = 0.0;
    }
    // Compute element coordinates 
    int nP = eptr[shearElem+1]-eptr[shearElem];
    for (int i = eptr[shearElem]; i < eptr[shearElem+1]; ++i) {
      maxLocationAndTime[0] += coordinates[connectivity[i]*ndim];
      maxLocationAndTime[1] += coordinates[connectivity[i]*ndim+1];
      maxLocationAndTime[2] += coordinates[connectivity[i]*ndim+2];
    }
    for (int i = 0; i < ndim; ++i) {
      maxLocationAndTime[i] = maxLocationAndTime[i]/((double)nP);
    }
    maxLocationAndTime[3] = maxShearT;
    MPI_Send(maxLocationAndTime, 4, MPI_DOUBLE, 0, 7299, MPI_COMM_WORLD);
    MPI_Send(&global_eid[shearElem], 1, MPI_INT, 0, 7299, MPI_COMM_WORLD);
  }
  if (world_rank == 0) {
    double minRecv[4];
    double maxRecv[4];
    double shearRecv[4];
    int globalIDMax[3];
    MPI_Recv(minRecv, 4, MPI_DOUBLE, parStructMin.rank, 7297, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&globalIDMax[0], 1, MPI_INT, parStructMin.rank, 7297, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(maxRecv, 4, MPI_DOUBLE, parStructMax.rank, 7298, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&globalIDMax[1], 1, MPI_INT, parStructMax.rank, 7298, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(shearRecv, 4, MPI_DOUBLE, parStructSMax.rank, 7299, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&globalIDMax[2], 1, MPI_INT, parStructSMax.rank, 7299, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    Json::Value output;
    Json::Value vec(Json::arrayValue);

    vec.append(Json::Value(maxRecv[0]));
    vec.append(Json::Value(maxRecv[1]));
    vec.append(Json::Value(maxRecv[2]));
    output["principal-max-strain"]["value"] = maxStrain;
    output["principal-max-strain"]["location"] = vec;
    output["principal-max-strain"]["time"] = maxRecv[3];
    output["principal-max-strain"]["global-element-id"] = globalIDMax[1];

    vec[0] = minRecv[0];
    vec[1] = minRecv[1];
    vec[2] = minRecv[2];
    output["principal-min-strain"]["value"] = minStrain;
    output["principal-min-strain"]["location"] = vec;
    output["principal-min-strain"]["time"] = minRecv[3];
    output["principal-min-strain"]["global-element-id"] = globalIDMax[0];

    vec[0] = shearRecv[0];
    vec[1] = shearRecv[1];
    vec[2] = shearRecv[2];
    output["maximum-shear-strain"]["value"] = maxShear;
    output["maximum-shear-strain"]["location"] = vec;
    output["maximum-shear-strain"]["time"] = shearRecv[3];
    output["maximum-shear-strain"]["global-element-id"] = globalIDMax[2];

    output["output-file"] = "coarse_brain.pvd";

    Json::StreamWriterBuilder builder;
    builder["commentStyle"] = "None";
    builder["indentation"] = "  ";
    std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
    std::ofstream oFile("output_"+ uid + ".json");
    writer -> write(output, &oFile);
  }
}

void CalculateInjuryCriterions(void) {
  double currentStrainMaxElem, currentStrainMinElem, currentShearMaxElem;
  double currentStrainMax = 0.0, currentStrainMin = 0.0, currentShearMax = 0.0;
  for (int i = 0; i < nelements; i++) {
    CalculateMaximumPrincipalStrain(i, &currentStrainMaxElem, &currentStrainMinElem);
    currentShearMaxElem = currentStrainMaxElem-currentStrainMinElem;
    if (currentStrainMax < currentStrainMaxElem) {
      currentStrainMax = currentStrainMaxElem;
      maxElem = i;
    }
    if (currentStrainMin > currentStrainMinElem) {
      currentStrainMin = currentStrainMinElem;
      minElem = i;
    }
    if (currentShearMax < currentShearMaxElem) {
      currentShearMax = currentShearMaxElem;
      shearElem = i;
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
  if (currentShearMax > maxShear) {
    maxShearT = Time;
    maxShear = currentShearMax;
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

  parStructSMax.value = maxShear;
  parStructSMax.rank = world_rank;
  MPI_Allreduce(MPI_IN_PLACE, &parStructSMax, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
  maxShear = parStructSMax.value;
}
