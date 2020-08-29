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
void InitCustomPlot(const Json::Value& jsonInput);
void InitBoundaryCondition(const Json::Value& jsonInput);
void InitInjuryCriterion(void);
void updateBoundaryNeighbour(void);
void ApplyAccBoundaryConditions();
int getImpactID(std::string location);
void WriteOutputFile();
void CalculateInjuryCriterions(void);
void TransformMesh(const Json::Value& jsonInput);

/* Global Variables/Parameters */
double Time, dt;
int nSteps;
bool ImplicitStatic = false;
bool ImplicitDynamic = false;
bool ExplicitDynamic = true;
double ExplicitTimeStepReduction = 0.8;
double FailureTimeStep = 5e-8; // Set for max runtime of around 5 hrs on aws
int nPlotSteps = 50;
int nWriteSteps = 2000;

/* Global variables used only in this file */
int nodeIDtoPlot;
bool rankForCustomPlot;
/* Global variables for bc */
int *boundaryID = NULL;
int boundarySize;
double peakTime, tMax;
bool writeField = true;
/* Global variables for output */
int *outputNodeList;
int outputNodeCount = 0;
int *outputElemList;
int outputElemCount = 0;
MPI_Comm output_comm;
MPI_File outputFilePtr;
int output_size, output_rank;
double *outputNodeRigidDisp;
double *outputElemStress;

/* Variables to compute maximim and minimum strain */
double maxStrain = 0.0, minStrain = 0.0, maxShear = 0.0;
int maxElem = 0, minElem = 0, shearElem = 0;
double maxT = 0.0, minT = 0.0, maxShearT = 0.0;
struct {
  double value;
  int   rank;
} parStructMax, parStructMin, parStructSMax, parStructPSxSR; 
/* Variables for other injury metrics */
int *MPSgt15, *MPSgt30; // Variable to store if MPS exceeds 15 and 30
/* Variables for maximum Principal Strain times Strain Rate */
double *PS_Old, maxPSxSR, maxTimePSxSR;
int maxElemPSxSR;
int nElementsInjury, *elementIDInjury;
// part ID to exclude from injury computation
// 0 : Skull, 1 : CSF
const int injuryExcludePIDCount = 2;
const int injuryExcludePID[injuryExcludePIDCount] = {0, 1};

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
  Json::Value inputJson = InitFemTech(argc, argv);
  Json::Value simulationJson = inputJson["simulation"];

  std::string meshFile = simulationJson["mesh"].asString();
  tMax = simulationJson["maximum-time"].asDouble()/1000.0; // Convert to seconds
  if (!simulationJson["write-vtu"].empty()) {
    writeField = simulationJson["write-vtu"].asBool();
  }
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
  TransformMesh(simulationJson);
  InitCustomPlot(simulationJson);
  InitInjuryCriterion();

  // Initial settings for BC evaluations
  // Used if initial velocity and acceleration BC is to be set.
  InitBoundaryCondition(simulationJson);

  /* Write inital, undeformed configuration*/
  Time = 0.0;
  int plot_counter = 0;
  if (writeField) {
    WriteVTU(outputFileName.c_str(), plot_counter);
    WritePVD(outputFileName.c_str(), plot_counter);
  }
  stepTime[plot_counter] = Time;

  int time_step_counter = 0;
  /** Central Difference Method - Beta and Gamma */
  // double beta = 0;
  // double gamma = 0.5;

  ShapeFunctions();
  /*  Step-1: Calculate the mass matrix similar to that of belytschko. */
  AssembleLumpedMass();
  // Needs to be after shapefunctions
  CustomPlot();

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
  if (nsteps_plot == 0) {
    nsteps_plot = nSteps;
  }
  int nsteps_write = (int)(nSteps / nWriteSteps);
  if (nsteps_write == 0) {
    nsteps_write = nSteps;
  }

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
    int writeFileFlag = time_step_counter%nsteps_write;
    CheckEnergy(Time, writeFlag);

    CalculateInjuryCriterions();
    if (writeFileFlag == 0) {
      CustomPlot();
    }
    if (writeFlag == 0) {
      FILE_LOG_MASTER(INFO, "Time : %15.6e, dt=%15.6e, tmax : %15.6e", Time, dt, tMax);
      plot_counter = plot_counter + 1;

      if (plot_counter < MAXPLOTSTEPS) {
        stepTime[plot_counter] = Time;
        if (writeField) {
          FILE_LOG(INFO, "------ Plot %d: WriteVTU", plot_counter);
          WriteVTU(outputFileName.c_str(), plot_counter);
          WritePVD(outputFileName.c_str(), plot_counter);
        }
      }
      FILE_LOGMatrixRM(DEBUGLOG, displacements, nNodes, ndim, "Displacement Solution");
    }
    time_step_counter = time_step_counter + 1;
    dt = ExplicitTimeStepReduction * StableTimeStep();

    // Barrier not a must
    MPI_Barrier(MPI_COMM_WORLD);
  } // end explcit while loop
  // Write output if last step results not written
  int writeFileFlag = (time_step_counter-1)%nsteps_write;
  if (writeFileFlag != 0) {
    CustomPlot();
  }
  FILE_LOG_MASTER(INFO, "End of Iterative Loop");
  FILE_LOGMatrixRM(DEBUGLOG, displacements, nNodes, ndim, "Final Displacement Solution");

  WriteOutputFile();
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
  free1DArray(outputNodeList);
  free1DArray(outputElemList);
  free1DArray(outputNodeRigidDisp);
  // Free variables used for injury criterions
  free1DArray(MPSgt15);
  free1DArray(MPSgt30);
  free1DArray(PS_Old);
  free1DArray(elementIDInjury);

  if (rankForCustomPlot) {
    MPI_File_close(&outputFilePtr);
  }
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
      locV[j] = coordinates[index + j];
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
  // Rotate points to plot to obtain their rigid displacements
  if (rankForCustomPlot) {
    for (int i = 0; i < outputNodeCount; i++) {
      int index = outputNodeList[i] * ndim;
      int index1 = i * ndim;
      for (int j = 0; j < ndim; ++j) {
        locV[j] = coordinates[index + j];
      }
      V[0] = 0.0; V[1] = locV[0]; V[2] = locV[1]; V[3] = locV[2];
      quaternionRotate(V, R, Rinv, Vp); // Vp = RVR^{-1}
      for (int j = 0; j < ndim; ++j) {
        // Displacement
        outputNodeRigidDisp[index1+j] = Vp[j+1] - locV[j] + yInt[9+j];
      }
    }
  }
  return;
}

void InitCustomPlot(const Json::Value& jsonInput) {
  rankForCustomPlot = false;
  if (!jsonInput["output-nodes"].empty()) {
    int outputSize = jsonInput["output-nodes"].size();
    int *outputNodes = (int*)malloc(outputSize*sizeof(int));
    int *outputLocalNodes = (int*)malloc(outputSize*sizeof(int));
    int *containsNode = (int*)calloc(outputSize, sizeof(int));
    int *containsNodeCum = (int*)calloc(outputSize, sizeof(int));
    // Store input node list to array
    jsonToArrayInt(outputNodes, jsonInput["output-nodes"]);

    // Check if node is present in the current process
    for (int i = 0; i < outputSize; ++i) {
      int nodeItem = outputNodes[i];
      for (int j = 0; j < nNodes; ++j) {
        if (globalNodeID[j] == nodeItem) {
          containsNode[i] = 1;
          outputLocalNodes[i] = j;
          break;
        }
      }
    }
    // If multiple procs have same node use the lowest rank
    MPI_Scan(containsNode, containsNodeCum, outputSize, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    for (int i = 0; i < outputSize; ++i) {
      if (containsNode[i] == 1 && containsNodeCum[i] == 1) {
        outputNodeCount = outputNodeCount + 1;
      }
    }
    outputNodeList = (int*)malloc(outputNodeCount*sizeof(int));
    int currentIndex = 0;
    for (int i = 0; i < outputSize; ++i) {
      if (containsNode[i] == 1 && containsNodeCum[i] == 1) {
        outputNodeList[currentIndex] = outputLocalNodes[i];
        currentIndex = currentIndex + 1;
      }
    }
    assert(currentIndex == outputNodeCount);
    if (outputNodeCount > 0) {
      rankForCustomPlot = true;
    }
    free(outputNodes);
    free(outputLocalNodes);
    free(containsNode);
    free(containsNodeCum);
    for (int i = 0; i < outputNodeCount; ++i) {
      int nodeL = outputNodeList[i];
      FILE_LOG_SINGLE(INFO, "Node ID : %d (%15.9e, %15.9e, %15.9e)", globalNodeID[nodeL], coordinates[nodeL*ndim], coordinates[nodeL*ndim+1], coordinates[nodeL*ndim+2]);
    }
  }   
  if (!jsonInput["output-elements"].empty()) {
    int outputSize = jsonInput["output-elements"].size();
    int *outputElements = (int*)malloc(outputSize*sizeof(int));
    int *outputLocalElements = (int*)malloc(outputSize*sizeof(int));
    int *containsElement = (int*)calloc(outputSize, sizeof(int));
    // Store input node list to array
    jsonToArrayInt(outputElements, jsonInput["output-elements"]);

    // Check if element is present in the current process
    // TODO : Used the fact that global element ID is sorted and use binary
    // search
    for (int i = 0; i < outputSize; ++i) {
      int elemItem = outputElements[i];
      for (int j = 0; j < nelements; ++j) {
        if (global_eid[j] == elemItem) {
          containsElement[i] = 1;
          outputLocalElements[i] = j;
          outputElemCount = outputElemCount + 1;
          break;
        }
      }
    }
    if (outputElemCount > 0) {
      rankForCustomPlot = true;
      outputElemList = (int*)malloc(sizeof(int)*outputElemCount);
      int currentIndex = 0;
      for (int i = 0; i < outputSize; ++i) {
        if (containsElement[i]) {
          outputElemList[currentIndex] = outputLocalElements[i];
          currentIndex = currentIndex + 1;
        }
      }
      assert(currentIndex == outputElemCount);
    }
    free(outputElements);
    free(outputLocalElements);
    free(containsElement);
  }
  // Create communicator with process with output nodes
  if (rankForCustomPlot) {
    int color = 1;
    MPI_Comm_split(MPI_COMM_WORLD, color, world_rank, &output_comm);
  } else {
    MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, world_rank, &output_comm);
  }
  if (rankForCustomPlot) {
    // Get the number of processes in output communicator
    MPI_Comm_size(output_comm, &output_size);
    // Get the rank of the process in output communicator
    MPI_Comm_rank(output_comm, &output_rank);
    // Open output file and write header
    std::string outFileStr = "plot_"+ uid + ".dat";
    MPI_Info infoin;
    MPI_Info_create(&infoin);
    MPI_Info_set(infoin, "access_style", "write_once,random");
    const char *outFileName = outFileStr.c_str();

    int err;
    err = MPI_File_open(output_comm, outFileName, MPI_MODE_EXCL|MPI_MODE_WRONLY|MPI_MODE_CREATE, infoin, &outputFilePtr);
    if (err != MPI_SUCCESS) {
      if (output_rank == 0) {
        MPI_File_delete(outFileName, MPI_INFO_NULL);
      }
      MPI_Barrier(output_comm);
      err = MPI_File_open(output_comm, outFileName, MPI_MODE_EXCL|MPI_MODE_WRONLY|MPI_MODE_CREATE, infoin, &outputFilePtr);
      if (err != MPI_SUCCESS) {
        FILE_LOG_SINGLE(ERROR, "Unable to open file to write co-ordinates");
        TerminateFemTech(10);
      }
    }
    // Write the header 
    char *outputHeader = (char*)malloc(sizeof(char)*(10+outputNodeCount*60+outputElemCount*40));
    outputHeader[0] = 0;
    if (output_rank == 0) {
      sprintf(outputHeader, "#Time");
    }       
    for (int i = 0; i < outputNodeCount; ++i) {
      int globalN = globalNodeID[outputNodeList[i]];
      sprintf(outputHeader, "%s Node%08d-DispX Node%08d-DispY Node%08d-DispZ", outputHeader, globalN, globalN, globalN);
    }
    for (int i = 0; i < outputElemCount; ++i) {
      int globalE = global_eid[outputElemList[i]];
      sprintf(outputHeader, "%s Elem%08d-StreP Elem%08d-StreS", outputHeader, globalE, globalE);
    }
    if (output_rank == (output_size-1)) {
      sprintf(outputHeader, "%s\n", outputHeader);
    }
    MPI_File_write_ordered(outputFilePtr, outputHeader, strlen(outputHeader), MPI_CHAR, MPI_STATUS_IGNORE);
    free(outputHeader);
    outputNodeRigidDisp = (double*)calloc(outputNodeCount*ndim, sizeof(double));
  }
  return;
}

void CustomPlot() {
  if (rankForCustomPlot) {
    // Write the output 
    char *output = (char*)malloc(sizeof(char)*(17+outputNodeCount*51+outputElemCount*34));
    output[0] = 0;
    if (output_rank == 0) {
      sprintf(output, "%15.9e", Time);
    }       
    for (int i = 0; i < outputNodeCount; ++i) {
      unsigned int plotNode = outputNodeList[i]*ndim;
      unsigned int plotIndex = i*ndim;
      double xCord = displacements[plotNode]-outputNodeRigidDisp[plotIndex];
      double yCord = displacements[plotNode+1]-outputNodeRigidDisp[plotIndex+1];
      double zCord = displacements[plotNode+2]-outputNodeRigidDisp[plotIndex+2];
      sprintf(output, "%s %15.9e %15.9e %15.9e", output, xCord, yCord, zCord);
    }
    double currentStrainMaxElem, currentStrainMinElem, currentShearMaxElem;
    for (int i = 0; i < outputElemCount; ++i) {
      unsigned int plotElem = outputElemList[i];
      // Recalculation to decouple from Injury criterion
      if (materialID[pid[plotElem]] != 0) {
        CalculateMaximumPrincipalStrain(plotElem, &currentStrainMaxElem, \
            &currentStrainMinElem, &currentShearMaxElem);
      } else {
        currentStrainMaxElem = 0.0;
        currentShearMaxElem = 0.0;
      }
      sprintf(output, "%s %15.9e %15.9e", output, \
          currentStrainMaxElem, currentShearMaxElem);
    }
    if (output_rank == (output_size-1)) {
      sprintf(output, "%s\n", output);
    }
    MPI_File_write_ordered(outputFilePtr, output, strlen(output), MPI_CHAR, MPI_STATUS_IGNORE);
    free(output);
  }
  return;
}

void InitBoundaryCondition(const Json::Value& jsonInput) {
  static const double gC = 9.81;
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
      // linAccXv[i] = gC*linAccXv[i];
      linAccXt[i] = 0.001*linAccXt[i];
    }
    for (int i = 0; i < linAccYSize; ++i) {
      // linAccYv[i] = gC*linAccYv[i];
      linAccYt[i] = 0.001*linAccYt[i];
    }
    for (int i = 0; i < linAccZSize; ++i) {
      // linAccZv[i] = gC*linAccZv[i];
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
    peakTime = jsonInput["time-peak-acceleration"].asDouble()/1000.0; //Convert to second
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
        double ds = impactNodeCoord[i];
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
  // Convert angular to linear frame
  if (!jsonInput["angular-to-linear-frame"].empty()) {
    double factor[ndim];
    int index[ndim];
    for (int i = 0; i < ndim; ++i) {
      index[i] = -10;
      std::string tr = jsonInput["angular-to-linear-frame"][i].asString();
      if (tr.length() == 2) {
        switch(tr.at(1)) {
          case 'x' : index[i] = 0;
                     break;
          case 'y' : index[i] = 1;
                     break;
          case 'z' : index[i] = 2;
                     break;
          default : FILE_LOG_MASTER(ERROR, "Transformation has to be x/y/z");
                    TerminateFemTech(12);
                    break;
        }
        switch(tr.at(0)) {
          case '-' : factor[i] = -1;
                     break;
          case '+' : factor[i] = 1;
                     break;
          default : FILE_LOG_MASTER(ERROR, "Prefix of axis not -/+, check input transformation");
                    TerminateFemTech(12);
                    break;
        }
      } else {
        if (tr.length() == 1) {
          factor[i] = 1;
          switch(tr.at(0)) {
            case 'x' : index[i] = 0;
                      break;
            case 'y' : index[i] = 1;
                      break;
            case 'z' : index[i] = 2;
                      break;
            default : FILE_LOG_MASTER(ERROR, "Transformation has to be x/y/z");
                      TerminateFemTech(12);
                      break;
          }
        } else {
          FILE_LOG_MASTER(ERROR, "Error in angular to linear frame transformation. Please check input file");
          TerminateFemTech(12);
        }
      }
    }
    // Validate the input
    int indexSum = index[0] + index[1] + index[2];
    if (indexSum != 3) {
      FILE_LOG_MASTER(ERROR, "Error in angular to linear frame transformation. x, y and z not present");
      TerminateFemTech(12);
    }
    if ((index[0] != 0) && (index[1] != 0) && (index[2] != 0)) {
      FILE_LOG_MASTER(ERROR, "Error in angular to linear frame transformation. x not present");
      TerminateFemTech(12);
    }
    // Transform all angular acceleration values
    // Transforming individual values rather than pointer rotation for
    // readability
    double accSize[ndim];
    accSize[index[0]] = angAccXSize;
    accSize[index[1]] = angAccYSize;
    accSize[index[2]] = angAccZSize;
    double *angAccTNew[ndim], *angAccVNew[ndim];
    angAccTNew[0] = (double*)malloc(accSize[0]*sizeof(double));
    angAccVNew[0] = (double*)malloc(accSize[0]*sizeof(double));
    angAccTNew[1] = (double*)malloc(accSize[1]*sizeof(double));
    angAccVNew[1] = (double*)malloc(accSize[1]*sizeof(double));
    angAccTNew[2] = (double*)malloc(accSize[2]*sizeof(double));
    angAccVNew[2] = (double*)malloc(accSize[2]*sizeof(double));

    int transformedIndex = index[0];
    for (int i = 0; i < angAccXSize; ++i) {
      angAccTNew[transformedIndex][i] = angAccXt[i];
      angAccVNew[transformedIndex][i] = angAccXv[i]*factor[0];
    }
    transformedIndex = index[1];
    for (int i = 0; i < angAccYSize; ++i) {
      angAccTNew[transformedIndex][i] = angAccYt[i];
      angAccVNew[transformedIndex][i] = angAccYv[i]*factor[1];
    }
    transformedIndex = index[2];
    for (int i = 0; i < angAccZSize; ++i) {
      angAccTNew[transformedIndex][i] = angAccZt[i];
      angAccVNew[transformedIndex][i] = angAccZv[i]*factor[2];
    }

    angAccXSize = accSize[0];
    angAccYSize = accSize[1];
    angAccZSize = accSize[2];
    free(angAccXt);
    free(angAccXv);
    free(angAccYt);
    free(angAccYv);
    free(angAccZt);
    free(angAccZv);
    angAccXt = angAccTNew[0];
    angAccXv = angAccVNew[0];
    angAccYt = angAccTNew[1];
    angAccYv = angAccVNew[1];
    angAccZt = angAccTNew[2];
    angAccZv = angAccVNew[2];
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
  const std::map<std::string, int> pointMap{ {"top-right", 5285}, \
    {"top-left", 1754}, {"front-low", 2576}, {"front-high", 1968}, \
    {"right-high", 5677}, {"bottom-front", 2576}, {"top-front", 1968}, \
    {"top-rear", 1874}};
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

void WriteOutputFile() {
  // Find the gloabl min and max strain
  parStructMax.value = maxStrain;
  parStructMax.rank = world_rank;
  MPI_Allreduce(MPI_IN_PLACE, &parStructMax, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

  parStructMin.value = minStrain;
  parStructMin.rank = world_rank;
  MPI_Allreduce(MPI_IN_PLACE, &parStructMin, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);

  parStructSMax.value = maxShear;
  parStructSMax.rank = world_rank;
  MPI_Allreduce(MPI_IN_PLACE, &parStructSMax, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

  // Find the gloabl max PSxSR
  parStructPSxSR.value = maxPSxSR;
  parStructPSxSR.rank = world_rank;
  MPI_Allreduce(MPI_IN_PLACE, &parStructPSxSR, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

  // Compute the volume of all parts
  double volumePart[nPIDglobal];
  double elementVolume[nelements]; // For MPS computation
  for (int i = 0; i < nPIDglobal; ++i) {
    volumePart[i] = 0.0;
  }
  computePartVolume(volumePart, elementVolume);

  // Compute Volume with MPS > 15, 30
  double MPSgt15Volume = 0.0, MPSgt30Volume = 0.0;
  for (int i = 0; i < nElementsInjury; ++i) {
    if (MPSgt15[i]) {
      MPSgt15Volume = MPSgt15Volume + elementVolume[elementIDInjury[i]];
      if (MPSgt30[i]) {
        MPSgt30Volume = MPSgt30Volume + elementVolume[elementIDInjury[i]];
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &MPSgt15Volume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &MPSgt30Volume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

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
  if (parStructPSxSR.rank == world_rank) {
    double maxLocationAndTime[4];
    for (int i = 0; i < 4; ++i) {
      maxLocationAndTime[i] = 0.0;
    }
    // Compute element coordinates 
    int nP = eptr[maxElemPSxSR+1]-eptr[maxElemPSxSR];
    for (int i = eptr[maxElemPSxSR]; i < eptr[maxElemPSxSR+1]; ++i) {
      maxLocationAndTime[0] += coordinates[connectivity[i]*ndim];
      maxLocationAndTime[1] += coordinates[connectivity[i]*ndim+1];
      maxLocationAndTime[2] += coordinates[connectivity[i]*ndim+2];
    }
    for (int i = 0; i < ndim; ++i) {
      maxLocationAndTime[i] = maxLocationAndTime[i]/((double)nP);
    }
    maxLocationAndTime[3] = maxTimePSxSR;
    MPI_Send(maxLocationAndTime, 4, MPI_DOUBLE, 0, 7300, MPI_COMM_WORLD);
    MPI_Send(&global_eid[maxElemPSxSR], 1, MPI_INT, 0, 7300, MPI_COMM_WORLD);
  }

  if (world_rank == 0) {
    double minRecv[4];
    double maxRecv[4];
    double shearRecv[4];
    double maxPSxSrRecv[4];
    int globalIDMax[4];
    MPI_Recv(minRecv, 4, MPI_DOUBLE, parStructMin.rank, 7297, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&globalIDMax[0], 1, MPI_INT, parStructMin.rank, 7297, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(maxRecv, 4, MPI_DOUBLE, parStructMax.rank, 7298, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&globalIDMax[1], 1, MPI_INT, parStructMax.rank, 7298, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(shearRecv, 4, MPI_DOUBLE, parStructSMax.rank, 7299, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&globalIDMax[2], 1, MPI_INT, parStructSMax.rank, 7299, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(maxPSxSrRecv, 4, MPI_DOUBLE, parStructPSxSR.rank, 7300, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&globalIDMax[3], 1, MPI_INT, parStructPSxSR.rank, 7300, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // Excluded PID hardcoded in CSDM-15 computation
    int excludePID[2] = {0, 1};
    double totalVolume = 0.0;
    for (int i = 0; i < nPIDglobal; ++i) {
      if (i == excludePID[0] || i == excludePID[1]) {
        continue;
      }
      totalVolume += volumePart[i];
    }

    Json::Value output;
    Json::Value vec(Json::arrayValue);

    vec.append(Json::Value(maxRecv[0]));
    vec.append(Json::Value(maxRecv[1]));
    vec.append(Json::Value(maxRecv[2]));
    output["principal-max-strain"]["value"] = parStructMax.value;
    output["principal-max-strain"]["location"] = vec;
    output["principal-max-strain"]["time"] = maxRecv[3];
    output["principal-max-strain"]["global-element-id"] = globalIDMax[1];

    vec[0] = minRecv[0];
    vec[1] = minRecv[1];
    vec[2] = minRecv[2];
    output["principal-min-strain"]["value"] = parStructMin.value;
    output["principal-min-strain"]["location"] = vec;
    output["principal-min-strain"]["time"] = minRecv[3];
    output["principal-min-strain"]["global-element-id"] = globalIDMax[0];

    vec[0] = shearRecv[0];
    vec[1] = shearRecv[1];
    vec[2] = shearRecv[2];
    output["maximum-shear-strain"]["value"] = parStructSMax.value;
    output["maximum-shear-strain"]["location"] = vec;
    output["maximum-shear-strain"]["time"] = shearRecv[3];
    output["maximum-shear-strain"]["global-element-id"] = globalIDMax[2];

    vec[0] = maxPSxSrRecv[0];
    vec[1] = maxPSxSrRecv[1];
    vec[2] = maxPSxSrRecv[2];
    output["maximum-PSxSR"]["value"] = parStructPSxSR.value;
    output["maximum-PSxSR"]["location"] = vec;
    output["maximum-PSxSR"]["time"] = maxPSxSrRecv[3];
    output["maximum-PSxSR"]["global-element-id"] = globalIDMax[3];

    output["output-file"] = "coarse_brain.pvd";

    // Write center of mass co-ordinates in JSON
    vec[0] = cm[0];
    vec[1] = cm[1];
    vec[2] = cm[2];
    output["center-of-mass"] = vec;

    // Write brain volume in output.json
    output["brain-volume"] = totalVolume;

    // Write CSDM-15 and 30
    output["csdm-15"] = MPSgt15Volume/totalVolume;
    output["csdm-30"] = MPSgt30Volume/totalVolume;

    Json::StreamWriterBuilder builder;
    builder["commentStyle"] = "None";
    builder["indentation"] = "  ";
    std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
    std::ofstream oFile(uid + "_output.json");
    writer -> write(output, &oFile);
  }
}

void InitInjuryCriterion(void) {
  // Allocate variables for injury metrics
  nElementsInjury = 0;

  for (int i = 0; i < nelements; i++) {
    bool include = true;
    int elementPID = pid[i];
    for (int j = 0; j < injuryExcludePIDCount; ++j) {
      if (elementPID == injuryExcludePID[j]) {
        include = false;
        break;
      }
    }
    if (include) {
      nElementsInjury += 1;
    }
  }
  elementIDInjury = (int*)malloc(nElementsInjury*sizeof(int));
  int count = 0;
  for (int i = 0; i < nelements; i++) {
    bool include = true;
    int elementPID = pid[i];
    for (int j = 0; j < injuryExcludePIDCount; ++j) {
      if (elementPID == injuryExcludePID[j]) {
        include = false;
        break;
      }
    }
    if (include) {
      elementIDInjury[count] = i;
      count += 1;
    }
  }

  MPSgt15 = (int*)calloc(nElementsInjury, sizeof(int));
  MPSgt30 = (int*)calloc(nElementsInjury, sizeof(int));
  PS_Old = (double*)malloc(nElementsInjury*sizeof(double));
}

void CalculateInjuryCriterions(void) {
  double currentStrainMaxElem, currentStrainMinElem, currentShearMaxElem;
  double PSR, PSxSR;
  for (int j = 0; j < nElementsInjury; j++) {
    int i = elementIDInjury[j];
    CalculateMaximumPrincipalStrain(i, &currentStrainMaxElem, \
        &currentStrainMinElem, &currentShearMaxElem);
    if (maxStrain < currentStrainMaxElem) {
      maxStrain = currentStrainMaxElem;
      maxElem = i;
      maxT = Time;
    }
    if (minStrain > currentStrainMinElem) {
      minStrain = currentStrainMinElem;
      minElem = i;
      minT = Time;
    }
    if (maxShear < currentShearMaxElem) {
      maxShear = currentShearMaxElem;
      shearElem = i;
      maxShearT = Time;
    }
    // Calculate MPS 15
    if (!MPSgt15[j]) {
      if (currentStrainMaxElem > 0.15) {
        MPSgt15[j] = 1;
      }
    }
    // Calculate MPS 30
    if (!MPSgt30[j]) {
      if (currentStrainMaxElem > 0.30) {
        MPSgt30[j] = 1;
      }
    }
    // Compute maxPSxSR
    // Compute principal strain rate using first order backaward distance
    PSR = (currentStrainMaxElem - PS_Old[j])/dt;
    // TODO(Anil) : Sould absolute value be used for PSR ?
    PSxSR = currentStrainMaxElem*PSR; 
    if (maxPSxSR < PSxSR) {
      maxPSxSR = PSxSR;
      maxElemPSxSR = i;
      maxTimePSxSR = Time;
    }
    PS_Old[j] = currentStrainMaxElem;
  } // For loop over elements included for injury
}

void TransformMesh(const Json::Value& jsonInput) {
  double cordAngularSensor[3];
  // Check if mesh-transformation key is present in the input json file
  // If present transform mesh, otherwise subtract center of mass
  if (!jsonInput["head-cg"].empty()) {
    cm[0] = jsonInput["head-cg"][0].asDouble();
    cm[1] = jsonInput["head-cg"][1].asDouble();
    cm[2] = jsonInput["head-cg"][2].asDouble();
    FILE_LOG_MASTER(INFO, "Center of Mass used : (%15.9e, %15.9e, %15.9e)", cm[0], cm[1], cm[2]);
  } else {
    GetBodyCenterofMass(cm);
    FILE_LOG_MASTER(INFO, "Brain+Skull Center of Mass used : (%15.9e, %15.9e, %15.9e)", cm[0], cm[1], cm[2]);
  }
  if (!jsonInput["angular-sensor-position"].empty()) {
    cordAngularSensor[0] = jsonInput["angular-sensor-position"][0].asDouble();
    cordAngularSensor[1] = jsonInput["angular-sensor-position"][1].asDouble();
    cordAngularSensor[2] = jsonInput["angular-sensor-position"][2].asDouble();
    FILE_LOG_MASTER(INFO, "Angular sensor data about : (%15.9e, %15.9e, %15.9e)", \
        cordAngularSensor[0], cordAngularSensor[1], cordAngularSensor[2]);
  } else {
    for (unsigned int i = 0; i < ndim; ++i) {
      cordAngularSensor[i] = cm[i];
    }
    FILE_LOG_MASTER(INFO, "Angular sensor data about center of mass");
  }
  if (!jsonInput["mesh-transformation"].empty()) {
    double factor[ndim];
    int index[ndim];
    for (int i = 0; i < ndim; ++i) {
      index[i] = -10;
      std::string tr = jsonInput["mesh-transformation"][i].asString();
      if (tr.length() == 2) {
        switch(tr.at(1)) {
          case 'x' : index[i] = 0;
                     break;
          case 'y' : index[i] = 1;
                     break;
          case 'z' : index[i] = 2;
                     break;
          default : FILE_LOG_MASTER(ERROR, "Mesh transformation has to be x/y/z");
                    TerminateFemTech(12);
                    break;
        }
        switch(tr.at(0)) {
          case '-' : factor[i] = -1;
                     break;
          case '+' : factor[i] = 1;
                     break;
          default : FILE_LOG_MASTER(ERROR, "Prefix of axis not -/+, check input mesh-transformation");
                    TerminateFemTech(12);
                    break;
        }
      } else {
        if (tr.length() == 1) {
          factor[i] = 1;
          switch(tr.at(0)) {
            case 'x' : index[i] = 0;
                      break;
            case 'y' : index[i] = 1;
                      break;
            case 'z' : index[i] = 2;
                      break;
            default : FILE_LOG_MASTER(ERROR, "Mesh transformation has to be x/y/z");
                      TerminateFemTech(12);
                      break;
          }
        } else {
          FILE_LOG_MASTER(ERROR, "Error in mesh transformation. Please check input file");
          TerminateFemTech(12);
        }
      }
    }
    // Validate the input
    int indexSum = index[0] + index[1] + index[2];
    if (indexSum != 3) {
      FILE_LOG_MASTER(ERROR, "Error in mesh transformation. x, y and z not present");
      TerminateFemTech(12);
    }
    if ((index[0] != 0) && (index[1] != 0) && (index[2] != 0)) {
      FILE_LOG_MASTER(ERROR, "Error in mesh transformation. x not present");
      TerminateFemTech(12);
    }
    // Transform center of mass
    double coordTemp[ndim];
    // Transform center of mass based on the transformation
    for (int i = 0; i < ndim; ++i) {
      int transformedIndex = index[i];
      coordTemp[transformedIndex] = cm[i]*factor[i];
    }
    // Assign cm array with transformed co-ordinates
    for (int i = 0; i < ndim; ++i) {
      cm[i] = coordTemp[i];
    }
    // Transform angular sensor position
    for (int i = 0; i < ndim; ++i) {
      int transformedIndex = index[i];
      coordTemp[transformedIndex] = cordAngularSensor[i]*factor[i];
    }
    for (int i = 0; i < ndim; ++i) {
      cordAngularSensor[i] = coordTemp[i];
    }
    // Transform all corodinates and subtract angular sensor position.
    // Transform all local co-ordinates
    for (int j = 0; j < nNodes; ++j) {
      for (int i = 0; i < ndim; ++i) {
        int transformedIndex = index[i];
        coordTemp[transformedIndex] = coordinates[i + ndim*j]*factor[i];
      }
      for (int i = 0; i < ndim; ++i) {
        coordinates[i + ndim*j] = coordTemp[i] - cordAngularSensor[i];
      }
    }
    // Get center of mass wrt angular sensor position
    for (int i = 0; i < ndim; ++i) {
      cm[i] = cm[i] - cordAngularSensor[i];
    }
  } else {
    // Subtract angular sensor position from co-ordinates
    for (int j = 0; j < nNodes; ++j) {
      for (int i = 0; i < ndim; ++i) {
        coordinates[i + ndim*j] = coordinates[i + ndim*j] - cordAngularSensor[i];
      }
    }
    // Get center of mass wrt angular sensor position
    for (int i = 0; i < ndim; ++i) {
      cm[i] = cm[i] - cordAngularSensor[i];
    }
  }
  // Ensure rest of the code is executed after the mesh transformation
  MPI_Barrier(MPI_COMM_WORLD);
}
