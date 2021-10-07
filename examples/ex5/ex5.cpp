#include "FemTech.h"
#include "jsonfuncs.h"
#include "utilities.h"
#include "gitbranch.h"

#include "json/writer.h"

#include <assert.h>
#include <fstream> // For JSON output
#include <map>
#include <math.h>
#include <memory> // For JSON output unique pointer
#include <string>

#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

/*Declare Functions*/
void CustomPlot();
void InitCustomPlot(const Json::Value &jsonInput);
double InitBoundaryCondition(const Json::Value &jsonInput);
void InitInjuryCriteria(void);
void updateBoundaryNeighbour(void);
void ApplyAccBoundaryConditions();
int getImpactID(std::string location);
void WriteOutputFile();
void CalculateInjuryCriteria(void);
void TransformMesh(const Json::Value &jsonInput);
void WriteMPS(void);

/* Global Variables/Parameters */
double Time = 0.0, dt;
int nSteps;
bool ImplicitStatic = false;
bool ImplicitDynamic = false;
bool ExplicitDynamic = true;

double dynamicDamping = 0.000;
double ExplicitTimeStepReduction = 0.8;
double FailureTimeStep = 1e-8; // Set for max runtime of around 5 hrs on aws
double MaxTimeStep = 1e-5;
int nPlotSteps = 50;
int nWriteSteps = 100;

/* Global variables used only in this file */
int nodeIDtoPlot;
bool rankForCustomPlot;
bool outputRelativeDisplacement = true;
/* Global variables for bc */
int *boundaryID = NULL;
int boundarySize;
double peakTime, tMax;
bool writeField = true;
bool computeInjuryFlag = true;
/* Global variables for output */
int *outputNodeList;
int outputNodeCount = 0;
int *outputElemList;
int outputElemCount = 0;
MPI_Comm output_comm;
MPI_File outputFilePtr;
int output_size, output_rank;
double *outputElemStress;

/* Variables for injury metrics */
/* Variables to compute maximum quatities over time and space */
/* 0. Maximum principal strain, 1. Maximum principal shear,
 * 2. Maximum principal strain times strain rate
 * 3. Minimum principal strain
 * */
const int maxQuantities = 4;
const std::string maxOutput[maxQuantities] = {
    "principal-max-strain", "maximum-shear-strain", "maximum-PSxSR",
    "principal-min-strain"};
const MPI_Op outputOperator[maxQuantities] = {MPI_MAXLOC, MPI_MAXLOC,
                                              MPI_MAXLOC, MPI_MINLOC};

double maxValue[maxQuantities], maxTime[maxQuantities];
int maxElem[maxQuantities];

/* Percentile variables */
/* 0. MPS - 95 Percentile
 * 1. MPSxSR - 95 Percentile
 */
const int percentileQuantities = 2;
int *percentileElements[percentileQuantities];
double percentileTime[percentileQuantities];
double percentileValue[percentileQuantities];
const std::string percentileTag[percentileQuantities] = {"MPS-95", "MPSxSR-95"};
/* Threshold variables */
/* 0. MPS > 3
 * 1. MPS > 5
 * 2. MPS > 10
 * 3. MPS > 15
 * 4. MPS > 30
 * 5. MPSR > 120 1/s
 * 6. MPSxSR > 28 1/s
 */
const int threshQuantities = 7;
int *thresholdElements[threshQuantities];
const std::string thresholdTag[threshQuantities] = {
    "CSDM-3", "CSDM-5", "CSDM-10", "CSDM-15", "CSDM-30", "MPSR-120", "MPSxSR-28"};

/* Variables for maximum Principal Strain times Strain Rate */
double *PS_Old = NULL;
double *PSxSRArray = NULL;
/* Store MPS, volume of each element */
double *elementMPS = NULL;
double *initialVolume = NULL;

int nElementsInjury = 0, *elementIDInjury = NULL;
// part ID to exclude from injury computation
// 0 : Skull, 1 : CSF
const int injuryExcludePIDCount = 2;
const int injuryExcludePID[injuryExcludePIDCount] = {0, 1};

// Variables to write injury criterion to Paraview output
const int outputCount = 8;
int *outputDataArray[outputCount];
const char *outputNames[outputCount] = {"CSDM-3", "CSDM-5",  "CSDM-10",  "CSDM-15",
                                        "CSDM-30", "MPSR-120", "MPSxSR-28",
                                        "MPS-95"};
const int outputDoubleCount = 1;
double *outputDoubleArray[outputDoubleCount];
const char *outputDoubleNames[outputDoubleCount] = {"MPS-95-Value"};

const int csdmCount = 5;
const double csdmLimits[csdmCount] = {0.03, 0.05, 0.1, 0.15, 0.3};

/* Variables used to store acceleration values */
int linAccXSize, linAccYSize, linAccZSize;
int angAccXSize, angAccYSize, angAccZSize;
int angVelXSize, angVelYSize, angVelZSize;
double *linAccXt, *linAccYt, *linAccZt;
double *linAccXv, *linAccYv, *linAccZv;
double *angAccXt, *angAccYt, *angAccZt;
double *angAccXv, *angAccYv, *angAccZv;
double *angVelXt, *angVelYt, *angVelZt;
double *angVelXv, *angVelYv, *angVelZv;
bool angularVelPrescribed = false;

/* Local co-ordinate system rotating with the skull for comparison with Hardy's
 * data. Use first 3 index to store x unit vector, next 3 for y unit vector  */
double unitVec[9];

/* Varibles for Local integrator */
const int nIntVar = 12;
typedef boost::array<double, nIntVar> state_type;
runge_kutta_dopri5<state_type> rk;
void computeDerivatives(const state_type &y, state_type &ydot, const double t);
state_type yInt, ydotInt;
double cm[3];
std::string outputFileName;

int main(int argc, char **argv) {
  // Initialize FemTech including logfile and MPI
  Json::Value inputJson = InitFemTech(argc, argv);
  Json::Value simulationJson = inputJson["simulation"];

  std::string meshFile = simulationJson["mesh"].asString();
  // Convert maximum time to seconds
  tMax = simulationJson["maximum-time"].asDouble() / 1000.0;
  if (!simulationJson["write-vtu"].empty()) {
    writeField = simulationJson["write-vtu"].asBool();
  }
  if (!simulationJson["compute-injury-criteria"].empty()) {
    computeInjuryFlag = simulationJson["compute-injury-criteria"].asBool();
  }
  FILE_LOG_MASTER(INFO, "Reading Mesh File : %s", meshFile.c_str());
  // Read Input Mesh file and equally partition elements among processes
  ReadInputFile(meshFile.c_str());
  size_t lastindex = meshFile.find_last_of(".");
  outputFileName = meshFile.substr(0, lastindex) + "_" + uid;
  // Read material properties before mesh partition to estimate
  // material type kernel compute intensity
  ReadMaterials();

  PartitionMesh();

  AllocateArrays();
  TransformMesh(simulationJson);
  InitCustomPlot(simulationJson);
  if (computeInjuryFlag) {
    InitInjuryCriteria();
  }

  // Initial settings for BC evaluations
  // Used if initial velocity and acceleration BC is to be set.
  Time = InitBoundaryCondition(simulationJson);

  /* Write inital, undeformed configuration*/
  int plot_counter = 0;
  if (writeField) {
    if (computeInjuryFlag) {
      WriteVTU(outputFileName.c_str(), plot_counter, outputDataArray, outputCount,
              outputNames, elementIDInjury, nElementsInjury, outputDoubleArray,
              outputDoubleCount, outputDoubleNames);
    } else {
      WriteVTU(outputFileName.c_str(), plot_counter);
    }
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
  // GetForce(); // Calculating the force term.

  /* Step-3: Calculate accelerations */
  // CalculateAccelerations();
  for (int i = 0; i < nDOF; i++) {
    accelerations[i] = 0.0;
  }

  nSteps = (int)((tMax - Time) / dt);
  int nsteps_plot = (int)(nSteps / nPlotSteps);
  if (nsteps_plot == 0) {
    nsteps_plot = nSteps;
  }
  int nsteps_write = (int)(nSteps / nWriteSteps);
  if (nsteps_write == 0) {
    nsteps_write = nSteps;
  }

  double t_n = 0.0;

  FILE_LOG_MASTER(INFO, "---------------------------------");
  FILE_LOG_MASTER(INFO, "Tmax : %15.6e, Initial dt : %15.6e", tMax, dt);
  FILE_LOG_MASTER(INFO, "nSteps = %d, nsteps_plot = %d", nSteps, nsteps_plot);
  FILE_LOG_MASTER(INFO, "-------------- Loop -------------");

  time_step_counter = time_step_counter + 1;


  /* Step-4: Time loop starts....*/
  while (Time < tMax) {
    t_n = Time;
    double t_np1 = Time + dt;
    double dtby2 = 0.5*dt;
    double t_nphalf = t_n + dtby2; // equ 6.2.1

    /* Step 6 Enforce boundary Conditions */
    ApplyAccBoundaryConditions();
    Time = t_np1; /*Update the time by adding full time step */
    /* Step 5 from Belytschko Box 6.1 - Update velocity */
    for (int i = 0; i < nDOF; i++) {
      if (!boundary[i]) {
        velocities_half[i] = velocities[i] + dtby2 * accelerations[i];
      }
    }

    // Store old displacements and accelerations for energy computation
    memcpy(displacements_prev, displacements, nDOF * sizeof(double));
    memcpy(accelerations_prev, accelerations, nDOF * sizeof(double));
    // Store internal external force from previous step to compute energy
    memcpy(fi_prev, fi, nDOF * sizeof(double));
    memcpy(fe_prev, fe, nDOF * sizeof(double));

    // update displacements for all nodes, including where velocity bc is set
    for (int i = 0; i < nDOF; i++) {
      displacements[i] = displacements[i] + dt * velocities_half[i];
    }

    /* Step - 8 from Belytschko Box 6.1 - Calculate net nodal force*/
    GetForce(); // Calculating the force term.
    /* Step - 9 from Belytschko Box 6.1 - Calculate Accelerations */
    CalculateAccelerations(); // Calculating the new accelerations from total
                              // nodal forces.

    /** Step- 10 - Second Partial Update of Nodal Velocities */
    for (int i = 0; i < nDOF; i++) {
      velocities[i] = velocities_half[i] + dtby2 * accelerations[i];
    }

    /** Step - 11 Checking* Energy Balance */
    int writeFlag = time_step_counter % nsteps_plot;
    int writeFileFlag = time_step_counter % nsteps_write;
    CheckEnergy(Time, writeFlag);

    if (computeInjuryFlag) {
      CalculateInjuryCriteria();
    }
    if (writeFileFlag == 0) {
      CustomPlot();
    }
    if (writeFlag == 0) {
      FILE_LOG_MASTER(INFO, "Time : %15.6e, dt=%15.6e, tmax : %15.6e", Time, dt,
                      tMax);
      plot_counter = plot_counter + 1;

      if (plot_counter < MAXPLOTSTEPS) {
        stepTime[plot_counter] = Time;
        if (writeField) {
          FILE_LOG(INFO, "------ Plot %d: WriteVTU", plot_counter);
          if (computeInjuryFlag) {
            WriteVTU(outputFileName.c_str(), plot_counter, outputDataArray, outputCount,
                    outputNames, elementIDInjury, nElementsInjury, outputDoubleArray,
                    outputDoubleCount, outputDoubleNames);
          } else {
            WriteVTU(outputFileName.c_str(), plot_counter);
          }
          WritePVD(outputFileName.c_str(), plot_counter);
        }
      }
    }
    time_step_counter = time_step_counter + 1;
    dt = ExplicitTimeStepReduction * StableTimeStep();

    // Barrier not a must
    MPI_Barrier(MPI_COMM_WORLD);
  } // end explcit while loop
  // Write output if last step results not written
  int writeFileFlag = (time_step_counter - 1) % nsteps_write;
  if (writeFileFlag != 0) {
    CustomPlot();
  }
  FILE_LOG_MASTER(INFO, "End of Iterative Loop");

  WriteOutputFile();
  if (computeInjuryFlag) {
    WriteMPS();
  }

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
  free1DArray(angVelXt);
  free1DArray(angVelYt);
  free1DArray(angVelZt);
  free1DArray(angVelXv);
  free1DArray(angVelYv);
  free1DArray(angVelZv);
  free1DArray(outputNodeList);
  free1DArray(outputElemList);
  // Free variables used for injury criteria
  free1DArray(PS_Old);
  free1DArray(elementIDInjury);
  free1DArray(PSxSRArray);
  for (int i = 0; i < threshQuantities; ++i) {
    free1DArray(thresholdElements[i]);
  }
  for (int i = 0; i < percentileQuantities; ++i) {
    free1DArray(percentileElements[i]);
  }
  free1DArray(elementMPS);
  free1DArray(initialVolume);

  if (rankForCustomPlot) {
    MPI_File_close(&outputFilePtr);
  }
  FinalizeFemTech();
  return 0;
}

void ApplyAccBoundaryConditions() {
  double r[4], R[4], Rinv[4], V[4], Vp[4]; // quaternions
  double omegaR[3], omega[3]; // vectors
  double locV[3];

  // Update to half time step to compute velocity
  rk.do_step(computeDerivatives, yInt, Time, 0.5*dt);
  r[0] = 0.0;
  r[1] = yInt[3];
  r[2] = yInt[4];
  r[3] = yInt[5];
  quaternionExp(r, R); // R = exp(r)
  quaternionInverse(R, Rinv);
  if (angularVelPrescribed) {
    double tLocal = Time + 0.5*dt;
    omega[0] = interpolateLinear(angVelXSize, angVelXt, angVelXv, tLocal);
    omega[1] = interpolateLinear(angVelYSize, angVelYt, angVelYv, tLocal);
    omega[2] = interpolateLinear(angVelZSize, angVelZt, angVelZv, tLocal);
  } else {
    omega[0] = yInt[0];
    omega[1] = yInt[1];
    omega[2] = yInt[2];
  }

  for (int i = 0; i < boundarySize; i++) {
    int index = boundaryID[i] * ndim;
    for (int j = 0; j < ndim; ++j) {
      locV[j] = coordinates[index + j];
    }
    V[0] = 0.0;
    V[1] = locV[0];
    V[2] = locV[1];
    V[3] = locV[2];
    quaternionRotate(V, R, Rinv, Vp); // Vp = RVR^{-1}
    crossProduct(omega, &(Vp[1]), omegaR);
    for (int j = 0; j < ndim; ++j) {
      velocities_half[index + j] = omegaR[j] + yInt[6 + j];
    }
  }
  // Move to full time step
  rk.do_step(computeDerivatives, yInt, Time+0.5*dt, 0.5*dt);
  // Rotate unit vectors to keep track of local co-ordinate system
  if (rankForCustomPlot) {
    r[0] = 0.0;
    r[1] = yInt[3];
    r[2] = yInt[4];
    r[3] = yInt[5];
    quaternionExp(r, R); // R = exp(r)
    quaternionInverse(R, Rinv);
    // Update unit vector
    for (int i = 0; i < 3; i++) {
      int index = i * ndim;
      V[0] = 0.0;
      V[1] = double(i==0);
      V[2] = double(i==1);
      V[3] = double(i==2);
      quaternionRotate(V, R, Rinv, Vp); // Vp = RVR^{-1}
      for (int j = 0; j < ndim; ++j) {
        // update unit vectors by rotating
        unitVec[index + j] = Vp[j + 1];
      }
    }
  }
  return;
}

void InitCustomPlot(const Json::Value &jsonInput) {
  rankForCustomPlot = false;
  if (!jsonInput["output-nodes"].empty()) {
    int outputSize = jsonInput["output-nodes"].size();
    int *outputNodes = (int *)malloc(outputSize * sizeof(int));
    int *outputLocalNodes = (int *)malloc(outputSize * sizeof(int));
    int *containsNode = (int *)calloc(outputSize, sizeof(int));
    int *containsNodeCum = (int *)calloc(outputSize, sizeof(int));
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
    MPI_Scan(containsNode, containsNodeCum, outputSize, MPI_INT, MPI_SUM,
             MPI_COMM_WORLD);
    for (int i = 0; i < outputSize; ++i) {
      if (containsNode[i] == 1 && containsNodeCum[i] == 1) {
        outputNodeCount = outputNodeCount + 1;
      }
    }
    outputNodeList = (int *)malloc(outputNodeCount * sizeof(int));
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
      FILE_LOG_SINGLE(INFO, "Node ID : %d (%15.9e, %15.9e, %15.9e)",
                      globalNodeID[nodeL], coordinates[nodeL * ndim],
                      coordinates[nodeL * ndim + 1],
                      coordinates[nodeL * ndim + 2]);
    }
  }
  if (!jsonInput["output-elements"].empty()) {
    int outputSize = jsonInput["output-elements"].size();
    int *outputElements = (int *)malloc(outputSize * sizeof(int));
    int *outputLocalElements = (int *)malloc(outputSize * sizeof(int));
    int *containsElement = (int *)calloc(outputSize, sizeof(int));
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
      outputElemList = (int *)malloc(sizeof(int) * outputElemCount);
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
    std::string outFileStr = "plot_" + uid + ".dat";
    MPI_Info infoin;
    MPI_Info_create(&infoin);
    MPI_Info_set(infoin, "access_style", "write_once,random");
    const char *outFileName = outFileStr.c_str();

    int err;
    err = MPI_File_open(output_comm, outFileName,
                        MPI_MODE_EXCL | MPI_MODE_WRONLY | MPI_MODE_CREATE,
                        infoin, &outputFilePtr);
    if (err != MPI_SUCCESS) {
      if (output_rank == 0) {
        MPI_File_delete(outFileName, MPI_INFO_NULL);
      }
      MPI_Barrier(output_comm);
      err = MPI_File_open(output_comm, outFileName,
                          MPI_MODE_EXCL | MPI_MODE_WRONLY | MPI_MODE_CREATE,
                          infoin, &outputFilePtr);
      if (err != MPI_SUCCESS) {
        FILE_LOG_SINGLE(ERROR, "Unable to open file to write co-ordinates");
        TerminateFemTech(10);
      }
    }
    // Write the header
    char *outputHeader = (char *)malloc(
        sizeof(char) * (10 + outputNodeCount * 60 + outputElemCount * 40));
    outputHeader[0] = 0;
    if (output_rank == 0) {
      sprintf(outputHeader, "#Time");
    }
    for (int i = 0; i < outputNodeCount; ++i) {
      int globalN = globalNodeID[outputNodeList[i]];
      sprintf(outputHeader, "%s Node%08d-DispX Node%08d-DispY Node%08d-DispZ",
              outputHeader, globalN, globalN, globalN);
    }
    for (int i = 0; i < outputElemCount; ++i) {
      int globalE = global_eid[outputElemList[i]];
      sprintf(outputHeader, "%s Elem%08d-StreP Elem%08d-StreS", outputHeader,
              globalE, globalE);
    }
    if (output_rank == (output_size - 1)) {
      sprintf(outputHeader, "%s\n", outputHeader);
    }
    MPI_File_write_ordered(outputFilePtr, outputHeader, strlen(outputHeader),
                           MPI_CHAR, MPI_STATUS_IGNORE);
    free(outputHeader);
    // Set the unit vector for local co-ordinates
    for (int i = 0; i < ndim2; ++i) {
      unitVec[i] = 0.0;
    }
    unitVec[0] = 1.0;
    unitVec[4] = 1.0;
    unitVec[8] = 1.0;
  }
  // Check if we should disable relative displacement
  if (!jsonInput["relative-displacement"].empty()) {
    outputRelativeDisplacement = jsonInput["relative-displacement"].asBool();
  }
  return;
}

void CustomPlot() {
  if (rankForCustomPlot) {
    // Write the output
    char *output = (char *)malloc(
        sizeof(char) * (17 + outputNodeCount * 51 + outputElemCount * 34));
    output[0] = 0;
    if (output_rank == 0) {
      sprintf(output, "%15.9e", Time);
    }
    for (int i = 0; i < outputNodeCount; ++i) {
      unsigned int plotNode = outputNodeList[i] * ndim;
      unsigned int plotIndex = i * ndim;
      double xCordLocal, yCordLocal, zCordLocal;
      if (outputRelativeDisplacement) {
        double xCord = displacements[plotNode] + coordinates[plotNode] - yInt[9];
        double yCord = displacements[plotNode + 1] + coordinates[plotNode + 1] - yInt[10];
        double zCord = displacements[plotNode + 2] + coordinates[plotNode + 2] - yInt[11];
        xCordLocal = xCord*unitVec[0] + yCord*unitVec[1] + zCord*unitVec[2] - coordinates[plotNode];
        yCordLocal = xCord*unitVec[3] + yCord*unitVec[4] + zCord*unitVec[5] - coordinates[plotNode + 1];
        zCordLocal = xCord*unitVec[6] + yCord*unitVec[7] + zCord*unitVec[8] - coordinates[plotNode + 2];
      } else {
        xCordLocal = displacements[plotNode];
        yCordLocal = displacements[plotNode + 1];
        zCordLocal = displacements[plotNode + 2];
      }
      sprintf(output, "%s %15.9e %15.9e %15.9e", output, xCordLocal, yCordLocal, zCordLocal);
    }
    double currentStrainMaxElem, currentStrainMinElem, currentShearMaxElem;
    for (int i = 0; i < outputElemCount; ++i) {
      unsigned int plotElem = outputElemList[i];
      // Recalculation to decouple from Injury criteria
      if (materialID[pid[plotElem]] != 0) {
        CalculateMaximumPrincipalStrain(plotElem, &currentStrainMaxElem,
                                        &currentStrainMinElem,
                                        &currentShearMaxElem);
      } else {
        currentStrainMaxElem = 0.0;
        currentShearMaxElem = 0.0;
      }
      sprintf(output, "%s %15.9e %15.9e", output, currentStrainMaxElem,
              currentShearMaxElem);
    }
    if (output_rank == (output_size - 1)) {
      sprintf(output, "%s\n", output);
    }
    MPI_File_write_ordered(outputFilePtr, output, strlen(output), MPI_CHAR,
                           MPI_STATUS_IGNORE);
    free(output);
  }
  return;
}

double InitBoundaryCondition(const Json::Value &jsonInput) {
  static const double gC = 9.81;
  double startTime = 0.0;
  // Read input JSON for acceleration values
  if (jsonInput["linear-acceleration"]
          .isObject()) { // Use time traces from file
    // If time array is common
    if (!jsonInput["time-all"].empty()) {
      unsigned int timeSize = 0;
      timeSize = jsonInput["time-all"].size();
      linAccXSize = jsonInput["linear-acceleration"]["xv"].size();
      assert(linAccXSize == timeSize);
      linAccYSize = jsonInput["linear-acceleration"]["yv"].size();
      assert(linAccYSize == timeSize);
      linAccZSize = jsonInput["linear-acceleration"]["zv"].size();
      assert(linAccZSize == timeSize);

      linAccXt = (double *)malloc(sizeof(double) * linAccXSize);
      linAccXv = (double *)malloc(sizeof(double) * linAccXSize);
      linAccYt = (double *)malloc(sizeof(double) * linAccYSize);
      linAccYv = (double *)malloc(sizeof(double) * linAccYSize);
      linAccZt = (double *)malloc(sizeof(double) * linAccZSize);
      linAccZv = (double *)malloc(sizeof(double) * linAccZSize);

      jsonToArray(linAccXt, jsonInput["time-all"]);
      jsonToArray(linAccXv, jsonInput["linear-acceleration"]["xv"]);
      jsonToArray(linAccYv, jsonInput["linear-acceleration"]["yv"]);
      jsonToArray(linAccZv, jsonInput["linear-acceleration"]["zv"]);
      memcpy(linAccYt, linAccXt, timeSize*sizeof(double));
      memcpy(linAccZt, linAccXt, timeSize*sizeof(double));

      if (!jsonInput["angular-velocity"].empty()) {
        angVelXSize = jsonInput["angular-velocity"]["xv"].size();
        assert(timeSize == angVelXSize);
        angVelYSize = jsonInput["angular-velocity"]["yv"].size();
        assert(timeSize == angVelYSize);
        angVelZSize = jsonInput["angular-velocity"]["zv"].size();
        assert(timeSize == angVelZSize);

        angVelXt = (double *)malloc(sizeof(double) * angVelXSize);
        angVelXv = (double *)malloc(sizeof(double) * angVelXSize);
        angVelYt = (double *)malloc(sizeof(double) * angVelYSize);
        angVelYv = (double *)malloc(sizeof(double) * angVelYSize);
        angVelZt = (double *)malloc(sizeof(double) * angVelZSize);
        angVelZv = (double *)malloc(sizeof(double) * angVelZSize);

        jsonToArray(angVelXv, jsonInput["angular-velocity"]["xv"]);
        jsonToArray(angVelYv, jsonInput["angular-velocity"]["yv"]);
        jsonToArray(angVelZv, jsonInput["angular-velocity"]["zv"]);
        memcpy(angVelXt, linAccXt, timeSize*sizeof(double));
        memcpy(angVelYt, linAccXt, timeSize*sizeof(double));
        memcpy(angVelZt, linAccXt, timeSize*sizeof(double));
        angularVelPrescribed = true;
      } else {
        if (jsonInput["angular-acceleration"].empty()) {
          // If both angular acceleration or angular velocity is not prescribed 
          // terminate the application
          FILE_LOG_SINGLE(ERROR, "Either angular acceleration or angular velocity has to be prescribed");
          TerminateFemTech(12);
        }
        angAccXSize = jsonInput["angular-acceleration"]["xv"].size();
        assert(timeSize == angAccXSize);
        angAccYSize = jsonInput["angular-acceleration"]["yv"].size();
        assert(timeSize == angAccYSize);
        angAccZSize = jsonInput["angular-acceleration"]["zv"].size();
        assert(timeSize == angAccZSize);

        angAccXt = (double *)malloc(sizeof(double) * angAccXSize);
        angAccXv = (double *)malloc(sizeof(double) * angAccXSize);
        angAccYt = (double *)malloc(sizeof(double) * angAccYSize);
        angAccYv = (double *)malloc(sizeof(double) * angAccYSize);
        angAccZt = (double *)malloc(sizeof(double) * angAccZSize);
        angAccZv = (double *)malloc(sizeof(double) * angAccZSize);

        jsonToArray(angAccXv, jsonInput["angular-acceleration"]["xv"]);
        jsonToArray(angAccYv, jsonInput["angular-acceleration"]["yv"]);
        jsonToArray(angAccZv, jsonInput["angular-acceleration"]["zv"]);
        memcpy(angAccXt, linAccXt, timeSize*sizeof(double));
        memcpy(angAccYt, linAccXt, timeSize*sizeof(double));
        memcpy(angAccZt, linAccXt, timeSize*sizeof(double));
      }
    } else {
      // Read linear acceleration and angular acceleration time traces
      linAccXSize = jsonInput["linear-acceleration"]["xt"].size();
      int tempSize = jsonInput["linear-acceleration"]["xv"].size();
      assert(tempSize == linAccXSize);
      linAccXt = (double *)malloc(sizeof(double) * linAccXSize);
      linAccXv = (double *)malloc(sizeof(double) * linAccXSize);
      linAccYSize = jsonInput["linear-acceleration"]["yt"].size();
      tempSize = jsonInput["linear-acceleration"]["yv"].size();
      assert(tempSize == linAccYSize);
      linAccYt = (double *)malloc(sizeof(double) * linAccYSize);
      linAccYv = (double *)malloc(sizeof(double) * linAccYSize);
      linAccZSize = jsonInput["linear-acceleration"]["zt"].size();
      tempSize = jsonInput["linear-acceleration"]["zv"].size();
      assert(tempSize == linAccYSize);
      linAccZt = (double *)malloc(sizeof(double) * linAccZSize);
      linAccZv = (double *)malloc(sizeof(double) * linAccZSize);
      jsonToArray(linAccXt, jsonInput["linear-acceleration"]["xt"]);
      jsonToArray(linAccXv, jsonInput["linear-acceleration"]["xv"]);
      jsonToArray(linAccYt, jsonInput["linear-acceleration"]["yt"]);
      jsonToArray(linAccYv, jsonInput["linear-acceleration"]["yv"]);
      jsonToArray(linAccZt, jsonInput["linear-acceleration"]["zt"]);
      jsonToArray(linAccZv, jsonInput["linear-acceleration"]["zv"]);

      if (!jsonInput["angular-velocity"].empty()) {
        angVelXSize = jsonInput["angular-velocity"]["xt"].size();
        tempSize = jsonInput["angular-velocity"]["xv"].size();
        assert(tempSize == angVelXSize);
        angVelXt = (double *)malloc(sizeof(double) * angVelXSize);
        angVelXv = (double *)malloc(sizeof(double) * angVelXSize);
        angVelYSize = jsonInput["angular-velocity"]["yt"].size();
        tempSize = jsonInput["angular-velocity"]["yv"].size();
        assert(tempSize == angVelYSize);
        angVelYt = (double *)malloc(sizeof(double) * angVelYSize);
        angVelYv = (double *)malloc(sizeof(double) * angVelYSize);
        angVelZSize = jsonInput["angular-velocity"]["zt"].size();
        tempSize = jsonInput["angular-velocity"]["zv"].size();
        assert(tempSize == angVelZSize);
        angVelZt = (double *)malloc(sizeof(double) * angVelZSize);
        angVelZv = (double *)malloc(sizeof(double) * angVelZSize);
        jsonToArray(angVelXt, jsonInput["angular-velocity"]["xt"]);
        jsonToArray(angVelXv, jsonInput["angular-velocity"]["xv"]);
        jsonToArray(angVelYt, jsonInput["angular-velocity"]["yt"]);
        jsonToArray(angVelYv, jsonInput["angular-velocity"]["yv"]);
        jsonToArray(angVelZt, jsonInput["angular-velocity"]["zt"]);
        jsonToArray(angVelZv, jsonInput["angular-velocity"]["zv"]);
        angularVelPrescribed = true;
      } else {
        if (jsonInput["angular-acceleration"].empty()) {
          // If both angular acceleration or angular velocity is not prescribed 
          // terminate the application
          FILE_LOG_SINGLE(ERROR, "Either angular acceleration or angular velocity has to be prescribed");
          TerminateFemTech(12);
        }
        angAccXSize = jsonInput["angular-acceleration"]["xt"].size();
        tempSize = jsonInput["angular-acceleration"]["xv"].size();
        assert(tempSize == angAccXSize);
        angAccXt = (double *)malloc(sizeof(double) * angAccXSize);
        angAccXv = (double *)malloc(sizeof(double) * angAccXSize);
        angAccYSize = jsonInput["angular-acceleration"]["yt"].size();
        tempSize = jsonInput["angular-acceleration"]["yv"].size();
        assert(tempSize == angAccYSize);
        angAccYt = (double *)malloc(sizeof(double) * angAccYSize);
        angAccYv = (double *)malloc(sizeof(double) * angAccYSize);
        angAccZSize = jsonInput["angular-acceleration"]["zt"].size();
        tempSize = jsonInput["angular-acceleration"]["zv"].size();
        assert(tempSize == angAccZSize);
        angAccZt = (double *)malloc(sizeof(double) * angAccZSize);
        angAccZv = (double *)malloc(sizeof(double) * angAccZSize);
        jsonToArray(angAccXt, jsonInput["angular-acceleration"]["xt"]);
        jsonToArray(angAccXv, jsonInput["angular-acceleration"]["xv"]);
        jsonToArray(angAccYt, jsonInput["angular-acceleration"]["yt"]);
        jsonToArray(angAccYv, jsonInput["angular-acceleration"]["yv"]);
        jsonToArray(angAccZt, jsonInput["angular-acceleration"]["zt"]);
        jsonToArray(angAccZv, jsonInput["angular-acceleration"]["zv"]);
      }
    }

    // Convert linear accelerations from g force to m/s^2
    // Convert time from milli-seconds to seconds
    for (int i = 0; i < linAccXSize; ++i) {
      // linAccXv[i] = gC*linAccXv[i];
      linAccXt[i] = 0.001 * linAccXt[i];
    }
    for (int i = 0; i < linAccYSize; ++i) {
      // linAccYv[i] = gC*linAccYv[i];
      linAccYt[i] = 0.001 * linAccYt[i];
    }
    for (int i = 0; i < linAccZSize; ++i) {
      // linAccZv[i] = gC*linAccZv[i];
      linAccZt[i] = 0.001 * linAccZt[i];
    }
    if (angularVelPrescribed) {
      for (int i = 0; i < angVelXSize; ++i) {
        angVelXt[i] = 0.001 * angVelXt[i];
      }
      for (int i = 0; i < angVelYSize; ++i) {
        angVelYt[i] = 0.001 * angVelYt[i];
      }
      for (int i = 0; i < angVelZSize; ++i) {
        angVelZt[i] = 0.001 * angVelZt[i];
      }
    } else {
      for (int i = 0; i < angAccXSize; ++i) {
        angAccXt[i] = 0.001 * angAccXt[i];
      }
      for (int i = 0; i < angAccYSize; ++i) {
        angAccYt[i] = 0.001 * angAccYt[i];
      }
      for (int i = 0; i < angAccZSize; ++i) {
        angAccZt[i] = 0.001 * angAccZt[i];
      }
    }
    startTime = linAccXt[0];
  } else { // Read maximum values from input file
    double accMax[3];
    accMax[0] = gC * jsonInput["linear-acceleration"][0].asDouble();
    accMax[1] = gC * jsonInput["linear-acceleration"][1].asDouble();
    accMax[2] = gC * jsonInput["linear-acceleration"][2].asDouble();
    peakTime = jsonInput["time-peak-acceleration"].asDouble() /
               1000.0; // Convert to second
    linAccXSize = 3;
    linAccYSize = 3;
    linAccZSize = 3;
    linAccXt = (double *)calloc(linAccXSize, sizeof(double));
    linAccXv = (double *)calloc(linAccXSize, sizeof(double));
    linAccYt = (double *)calloc(linAccYSize, sizeof(double));
    linAccYv = (double *)calloc(linAccYSize, sizeof(double));
    linAccZt = (double *)calloc(linAccZSize, sizeof(double));
    linAccZv = (double *)calloc(linAccZSize, sizeof(double));
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
      int nodeStatus = coordinateFromGlobalID(globalNodeID, impactPointID,
                                              nNodes, impactNodeCoord);
      // Find lowest rank process with node ID
      int *nodeWithID = (int *)malloc(world_size * sizeof(int));
      MPI_Allgather(&nodeStatus, 1, MPI_INT, nodeWithID, 1, MPI_INT,
                    MPI_COMM_WORLD);
      int nodeIDGlobal = -1;
      for (int i = 0; i < world_size; ++i) {
        if (nodeWithID[i]) {
          nodeIDGlobal = i;
          break;
        }
      }
      // Recieve node co-ordinates
      MPI_Bcast(impactNodeCoord, ndim, MPI_DOUBLE, nodeIDGlobal,
                MPI_COMM_WORLD);
      FILE_LOG_MASTER(INFO, "NodeID of impact : %d (%15.9e, %15.9e, %15.9e)",
                      impactPointID, impactNodeCoord[0], impactNodeCoord[1],
                      impactNodeCoord[2]);
      // Compute the axis of rotation
      double norm = 0.0;
      for (int i = 0; i < ndim; ++i) {
        double ds = impactNodeCoord[i];
        angNormal[i] = ds;
        norm += ds * ds;
      }
      norm = sqrt(norm);
      for (int i = 0; i < ndim; ++i) {
        angNormal[i] /= norm;
      }
      double angAccMax = jsonInput["angular-acceleration"].asDouble();
      accMax[0] = angAccMax * angNormal[0];
      accMax[1] = angAccMax * angNormal[1];
      accMax[2] = angAccMax * angNormal[2];
      free(nodeWithID);
    }
    angAccXSize = 3;
    angAccYSize = 3;
    angAccZSize = 3;
    angAccXt = (double *)calloc(angAccXSize, sizeof(double));
    angAccXv = (double *)calloc(angAccXSize, sizeof(double));
    angAccYt = (double *)calloc(angAccYSize, sizeof(double));
    angAccYv = (double *)calloc(angAccYSize, sizeof(double));
    angAccZt = (double *)calloc(angAccZSize, sizeof(double));
    angAccZv = (double *)calloc(angAccZSize, sizeof(double));
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
        switch (tr.at(1)) {
        case 'x':
          index[i] = 0;
          break;
        case 'y':
          index[i] = 1;
          break;
        case 'z':
          index[i] = 2;
          break;
        default:
          FILE_LOG_MASTER(ERROR, "Transformation has to be x/y/z");
          TerminateFemTech(12);
          break;
        }
        switch (tr.at(0)) {
        case '-':
          factor[i] = -1;
          break;
        case '+':
          factor[i] = 1;
          break;
        default:
          FILE_LOG_MASTER(ERROR,
                          "Prefix of axis not -/+, check input transformation");
          TerminateFemTech(12);
          break;
        }
      } else {
        if (tr.length() == 1) {
          factor[i] = 1;
          switch (tr.at(0)) {
          case 'x':
            index[i] = 0;
            break;
          case 'y':
            index[i] = 1;
            break;
          case 'z':
            index[i] = 2;
            break;
          default:
            FILE_LOG_MASTER(ERROR, "Transformation has to be x/y/z");
            TerminateFemTech(12);
            break;
          }
        } else {
          FILE_LOG_MASTER(ERROR, "Error in angular to linear frame "
                                 "transformation. Please check input file");
          TerminateFemTech(12);
        }
      }
    }
    // Validate the input
    int indexSum = index[0] + index[1] + index[2];
    if (indexSum != 3) {
      FILE_LOG_MASTER(ERROR, "Error in angular to linear frame transformation. "
                             "x, y and z not present");
      TerminateFemTech(12);
    }
    if ((index[0] != 0) && (index[1] != 0) && (index[2] != 0)) {
      FILE_LOG_MASTER(
          ERROR,
          "Error in angular to linear frame transformation. x not present");
      TerminateFemTech(12);
    }
    // Transform all angular acceleration/velocity values
    if (angularVelPrescribed) {
      // Transforming individual values rather than pointer rotation for
      // readability
      double velSize[ndim];
      velSize[index[0]] = angVelXSize;
      velSize[index[1]] = angVelYSize;
      velSize[index[2]] = angVelZSize;
      double *angVelTNew[ndim], *angVelVNew[ndim];
      angVelTNew[0] = (double *)malloc(velSize[0] * sizeof(double));
      angVelVNew[0] = (double *)malloc(velSize[0] * sizeof(double));
      angVelTNew[1] = (double *)malloc(velSize[1] * sizeof(double));
      angVelVNew[1] = (double *)malloc(velSize[1] * sizeof(double));
      angVelTNew[2] = (double *)malloc(velSize[2] * sizeof(double));
      angVelVNew[2] = (double *)malloc(velSize[2] * sizeof(double));

      int transformedIndex = index[0];
      for (int i = 0; i < angVelXSize; ++i) {
        angVelTNew[transformedIndex][i] = angVelXt[i];
        angVelVNew[transformedIndex][i] = angVelXv[i] * factor[0];
      }
      transformedIndex = index[1];
      for (int i = 0; i < angVelYSize; ++i) {
        angVelTNew[transformedIndex][i] = angVelYt[i];
        angVelVNew[transformedIndex][i] = angVelYv[i] * factor[1];
      }
      transformedIndex = index[2];
      for (int i = 0; i < angVelZSize; ++i) {
        angVelTNew[transformedIndex][i] = angVelZt[i];
        angVelVNew[transformedIndex][i] = angVelZv[i] * factor[2];
      }

      angVelXSize = velSize[0];
      angVelYSize = velSize[1];
      angVelZSize = velSize[2];
      free(angVelXt);
      free(angVelXv);
      free(angVelYt);
      free(angVelYv);
      free(angVelZt);
      free(angVelZv);
      angVelXt = angVelTNew[0];
      angVelXv = angVelVNew[0];
      angVelYt = angVelTNew[1];
      angVelYv = angVelVNew[1];
      angVelZt = angVelTNew[2];
      angVelZv = angVelVNew[2];
    } else {
      double accSize[ndim];
      // Transforming individual values rather than pointer rotation for
      // readability
      accSize[index[0]] = angAccXSize;
      accSize[index[1]] = angAccYSize;
      accSize[index[2]] = angAccZSize;
      double *angAccTNew[ndim], *angAccVNew[ndim];
      angAccTNew[0] = (double *)malloc(accSize[0] * sizeof(double));
      angAccVNew[0] = (double *)malloc(accSize[0] * sizeof(double));
      angAccTNew[1] = (double *)malloc(accSize[1] * sizeof(double));
      angAccVNew[1] = (double *)malloc(accSize[1] * sizeof(double));
      angAccTNew[2] = (double *)malloc(accSize[2] * sizeof(double));
      angAccVNew[2] = (double *)malloc(accSize[2] * sizeof(double));

      int transformedIndex = index[0];
      for (int i = 0; i < angAccXSize; ++i) {
        angAccTNew[transformedIndex][i] = angAccXt[i];
        angAccVNew[transformedIndex][i] = angAccXv[i] * factor[0];
      }
      transformedIndex = index[1];
      for (int i = 0; i < angAccYSize; ++i) {
        angAccTNew[transformedIndex][i] = angAccYt[i];
        angAccVNew[transformedIndex][i] = angAccYv[i] * factor[1];
      }
      transformedIndex = index[2];
      for (int i = 0; i < angAccZSize; ++i) {
        angAccTNew[transformedIndex][i] = angAccZt[i];
        angAccVNew[transformedIndex][i] = angAccZv[i] * factor[2];
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
      displacements[index + j] = 0.0;
      velocities[index + j] = 0.0;
      accelerations[index + j] = 0.0;
    }
  }
  return startTime;
}

void updateBoundaryNeighbour(void) {
  int *sendNodeBoundary;
  int *recvNodeBoundary;
  int totalNodeToSend = sendNeighbourCountCum[sendProcessCount];
  recvNodeBoundary = (int *)malloc(ndim * totalNodeToSend * sizeof(int));
  sendNodeBoundary = (int *)malloc(ndim * totalNodeToSend * sizeof(int));
  // Update array to send
  for (int i = 0; i < totalNodeToSend; ++i) {
    int nodeIndex = ndim * sendNodeIndex[i];
    sendNodeBoundary[ndim * i] = boundary[nodeIndex];
    sendNodeBoundary[ndim * i + 1] = boundary[nodeIndex + 1];
    sendNodeBoundary[ndim * i + 2] = boundary[nodeIndex + 2];
  }
  // Create send requests
  int boundaryTag = 2798;
  MPI_Request *requestListSend =
      (MPI_Request *)malloc(sizeof(MPI_Request) * sendProcessCount);
  for (int i = 0; i < sendProcessCount; ++i) {
    int process = sendProcessID[i];
    int location = sendNeighbourCountCum[i] * ndim;
    int size = ndim * sendNeighbourCount[i];
    MPI_Isend(&(sendNodeBoundary[location]), size, MPI_INT, process,
              boundaryTag, MPI_COMM_WORLD, &(requestListSend[i]));
  }
  // Create recv requests
  MPI_Request *requestListRecv =
      (MPI_Request *)malloc(sizeof(MPI_Request) * sendProcessCount);
  for (int i = 0; i < sendProcessCount; ++i) {
    int process = sendProcessID[i];
    int location = sendNeighbourCountCum[i] * ndim;
    int size = ndim * sendNeighbourCount[i];
    MPI_Irecv(&(recvNodeBoundary[location]), size, MPI_INT, process,
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
    int nodeIndex = ndim * sendNodeIndex[i];
    if (recvNodeBoundary[ndim * i]) {
      boundary[nodeIndex] = 1;
    }
    if (recvNodeBoundary[ndim * i + 1]) {
      boundary[nodeIndex + 1] = 1;
    }
    if (recvNodeBoundary[ndim * i + 2]) {
      boundary[nodeIndex + 2] = 1;
    }
  }
  free(requestListSend);
  free(requestListRecv);
  free(recvNodeBoundary);
  free(sendNodeBoundary);
}

int getImpactID(std::string location) {
  const std::map<std::string, int> pointMap{
      {"top-right", 5285},  {"top-left", 1754},   {"front-low", 2576},
      {"front-high", 1968}, {"right-high", 5677}, {"bottom-front", 2576},
      {"top-front", 1968},  {"top-rear", 1874}};
  if (pointMap.find(location) == pointMap.end()) {
    FILE_LOG_SINGLE(ERROR,
                    "Impact location not found, check value of impact-point");
    TerminateFemTech(11);
  }
  return pointMap.at(location);
}

void computeDerivatives(const state_type &y, state_type &ydot, const double t) {
  // ydot[0-2] : Store angular acceleration, alpha
  // ydot[3-5] : Store derivative of rotation quaternion generator, rdot
  // ydot[6-8] : Linear acceleration of the center of mass, acc
  // ydot[9-11] : Linear velocity of the center of mass, vel
  double omega[3];
  if (angularVelPrescribed) {
    // Computation of omega from alpha is not required as alpha is already
    // available
    ydot[0] = 0.0;
    ydot[1] = 0.0;
    ydot[2] = 0.0;
    // Interpolate omega
    omega[0] = interpolateLinear(angVelXSize, angVelXt, angVelXv, t);
    omega[1] = interpolateLinear(angVelYSize, angVelYt, angVelYv, t);
    omega[2] = interpolateLinear(angVelZSize, angVelZt, angVelZv, t);
  } else {
    // Storing alpha
    ydot[0] = interpolateLinear(angAccXSize, angAccXt, angAccXv, t);
    ydot[1] = interpolateLinear(angAccYSize, angAccYt, angAccYv, t);
    ydot[2] = interpolateLinear(angAccZSize, angAccZt, angAccZv, t);
    omega[0] = y[0];
    omega[1] = y[1];
    omega[2] = y[2];
  }
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
  r[0] = y[3];
  r[1] = y[4];
  r[2] = y[5];
  double *rdot = &(ydot[3]);
  double rMagnitude = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
  if (rMagnitude < 1e-10) {
    // rdot = 0.5*omega
    rdot[0] = 0.5 * omega[0];
    rdot[1] = 0.5 * omega[1];
    rdot[2] = 0.5 * omega[2];
  } else {
    double rCotR = rMagnitude / tan(rMagnitude);
    crossProduct(omega, r, rdot); // compute rdot = omega x r
    // Normalize r by its magnitude
    for (int i = 0; i < 3; ++i) {
      r[i] = r[i] / rMagnitude;
    }
    double rDotOmega = dotProduct3D(r, omega); // compure r.omega/||r||
    for (int i = 0; i < 3; ++i) {
      // Compute 0.5 omegaxr + omega*(0.5||r|| cot(||r||) + (1-||r||cot(||r||))
      // *(r.omega)*r/(2*||r||*||r||)
      rdot[i] =
          0.5 * (rdot[i] + rCotR * y[i] + (1.0 - rCotR) * rDotOmega * r[i]);
    }
  }
}

void WriteOutputFile() {
  Json::Value output;
  Json::Value vec(Json::arrayValue);

  // Write basic info that does not require any computation
  if (world_rank == 0) {
    output["output-file"] = outputFileName + ".pvd";

    vec.resize(ndim);
    // Write center of mass co-ordinates in JSON
    for (int i = 0; i < ndim; ++i) {
      vec[i] = cm[i];
    }
    output["center-of-mass"] = vec;
    output["code-branch"] = GIT_BRANCH;
    output["code-version"] = GIT_COMMIT_HASH;
  }
  if (computeInjuryFlag) {
    // Compute and write principal values
    double maxLocationAndTime[4];
    int globalElementID;
    struct {
      double value;
      int rank;
    } parStructMax;

    /* Calculate principal min and max strain location and send to master */
    // Find the gloabl max principal strain
    for (int i = 0; i < maxQuantities; ++i) {
      parStructMax.value = maxValue[i];
      parStructMax.rank = world_rank;
      MPI_Allreduce(MPI_IN_PLACE, &parStructMax, 1, MPI_DOUBLE_INT,
                    outputOperator[i], MPI_COMM_WORLD);
      if (parStructMax.rank == world_rank) {
        for (int j = 0; j < 4; ++j) {
          maxLocationAndTime[j] = 0.0;
        }
        int elementID = maxElem[i];
        // Compute element coordinates
        int nP = eptr[elementID + 1] - eptr[elementID];
        for (int j = eptr[elementID]; j < eptr[elementID + 1]; ++j) {
          maxLocationAndTime[0] += coordinates[connectivity[j] * ndim];
          maxLocationAndTime[1] += coordinates[connectivity[j] * ndim + 1];
          maxLocationAndTime[2] += coordinates[connectivity[j] * ndim + 2];
        }
        for (int j = 0; j < ndim; ++j) {
          maxLocationAndTime[j] = maxLocationAndTime[j] / ((double)nP);
        }
        maxLocationAndTime[3] = maxTime[i];
        globalElementID = global_eid[elementID];
        if (world_rank != 0) {
          MPI_Send(maxLocationAndTime, 4, MPI_DOUBLE, 0, 7297 + i,
                  MPI_COMM_WORLD);
          MPI_Send(&globalElementID, 1, MPI_INT, 0, 7297 + i, MPI_COMM_WORLD);
        }
      }
      if (world_rank == 0) {
        if (parStructMax.rank != 0) {
          MPI_Recv(maxLocationAndTime, 4, MPI_DOUBLE, parStructMax.rank, 7297 + i,
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          MPI_Recv(&globalElementID, 1, MPI_INT, parStructMax.rank, 7297 + i,
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        for (int j = 0; j < ndim; ++j) {
          vec[j] = maxLocationAndTime[j];
        }
        output[maxOutput[i]]["value"] = parStructMax.value;
        output[maxOutput[i]]["location"] = vec;
        output[maxOutput[i]]["time"] = maxLocationAndTime[3];
        output[maxOutput[i]]["global-element-id"] = globalElementID;
      }
    }

    // Compute the volume of all parts
    double volumePart[nPIDglobal];
    double elementVolume[nelements]; // For MPS computation
    for (int i = 0; i < nPIDglobal; ++i) {
      volumePart[i] = 0.0;
    }
    computePartVolume(volumePart, elementVolume);
    // Compute Volume of threshold quantities
    double thresholdVolume[threshQuantities];
    for (int i = 0; i < threshQuantities; ++i) {
      thresholdVolume[i] = 0.0;
    }
    for (int i = 0; i < nElementsInjury; ++i) {
      int e = elementIDInjury[i];
      double eV = elementVolume[e];
      // CSDM 3, 5, 10, 15, 30 volume computations
      for (unsigned int j = 0; j < csdmCount; ++j) {
        if (thresholdElements[j][i]) {
          thresholdVolume[j] += eV;
        } else {
          break;
        }
      }
      // MPSR > 120 volume
      if (thresholdElements[csdmCount][i]) {
        thresholdVolume[csdmCount] += eV;
      }
      // MPSxSR > 28 volume
      if (thresholdElements[csdmCount + 1][i]) {
        thresholdVolume[csdmCount + 1] += eV;
      }
    }
    // Calculate cumulative volume on all ranks
    MPI_Allreduce(MPI_IN_PLACE, thresholdVolume, 6, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    if (world_rank == 0) {
      // Excluded PID hardcoded for CSDM computation
      double totalVolume = 0.0;
      for (int i = 0; i < nPIDglobal; ++i) {
        bool include = true;
        for (int j = 0; j < injuryExcludePIDCount; ++j) {
          if (injuryExcludePID[j] == i) {
            include = false;
            break;
          }
        }
        if (include) {
          totalVolume += volumePart[i];
        }
      }
      // Write brain volume in output.json
      output["brain-volume"] = totalVolume;

      // Write threshold values
      for (int i = 0; i < threshQuantities; ++i) {
        output[thresholdTag[i]]["value"] = thresholdVolume[i] / totalVolume;
      }
      // Write percentile values
      for (int i = 0; i < percentileQuantities; ++i) {
        output[percentileTag[i]]["value"] = percentileValue[i];
        output[percentileTag[i]]["time"] = percentileTime[i];
      }
    }
    // Create list to write to output.json
    const int outputJsonCount = csdmCount + 1;
    std::string jsonOutputTag[outputJsonCount];
    int *outputJsonArray[outputJsonCount];
    for (int i = 0; i < csdmCount; ++i) {
      jsonOutputTag[i] = thresholdTag[i];
      outputJsonArray[i] = thresholdElements[i];
    }
    // Write MPS-95 elements
    jsonOutputTag[csdmCount] = percentileTag[0];
    outputJsonArray[csdmCount] = percentileElements[0];

    int *countPerProc, *cumCountPerProc, *fullElemList;
    if (world_rank == 0) {
      countPerProc = (int *)malloc(world_size * sizeof(int));
      cumCountPerProc = (int *)malloc(world_size * sizeof(int));
    }
    for (int i = 0; i < outputJsonCount; ++i) {
      int *list = outputJsonArray[i];
      int count = 0;
      for (int j = 0; j < nElementsInjury; ++j) {
        if (list[j]) {
          count = count + 1;
        }
      }
      int numElem = count;
      int *localElemList = (int *)malloc(numElem * sizeof(int));
      count = 0;
      if (numElem) {
        for (int j = 0; j < nElementsInjury; ++j) {
          if (list[j]) {
            int e = elementIDInjury[j];
            localElemList[count] = global_eid[e];
            count = count + 1;
          }
        }
      }
      // Get total count
      int totalElems = 0;
      MPI_Allreduce(&numElem, &totalElems, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      // Get recv count
      MPI_Gather(&numElem, 1, MPI_INT, countPerProc, 1, MPI_INT, 0,
                MPI_COMM_WORLD);

      if (world_rank == 0) {
        fullElemList = (int *)malloc(totalElems * sizeof(int));
        cumCountPerProc[0] = 0;
        for (int j = 1; j < world_size; ++j) {
          cumCountPerProc[j] = countPerProc[j - 1] + cumCountPerProc[j - 1];
        }
      }
      // Get full array
      if (totalElems) {
        MPI_Gatherv(localElemList, numElem, MPI_INT, fullElemList, countPerProc,
                    cumCountPerProc, MPI_INT, 0, MPI_COMM_WORLD);
      }
      if (world_rank == 0) {
        Json::Value elementList(Json::arrayValue);
        elementList.resize(totalElems);
        for (int j = 0; j < totalElems; ++j) {
          elementList[j] = fullElemList[j];
        }
        free(fullElemList);
        output[jsonOutputTag[i]]["global-element-id"] = elementList;
      }
      free(localElemList);
      MPI_Barrier(MPI_COMM_WORLD);
    }
    if (world_rank == 0) {
      free(countPerProc);
      free(cumCountPerProc);
    }
  }

  // Write output to file
  if (world_rank == 0) {
    Json::StreamWriterBuilder builder;
    builder["commentStyle"] = "None";
    builder["indentation"] = "  ";
    std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
    std::ofstream oFile(uid + "_output.json");
    writer->write(output, &oFile);
  }
  FILE_LOG_MASTER(INFO, "Json output file written");
}

void InitInjuryCriteria(void) {
  // Set all min and max quantities to zero
  for (int i = 0; i < maxQuantities; ++i) {
    maxValue[i] = 0.0;
    maxTime[i] = 0.0;
    maxElem[i] = 0;
  }

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
  elementIDInjury = (int *)malloc(nElementsInjury * sizeof(int));
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

  for (int i = 0; i < threshQuantities; ++i) {
    thresholdElements[i] = (int *)calloc(nElementsInjury, sizeof(int));
  }
  // Write percentile values
  for (int i = 0; i < percentileQuantities; ++i) {
    percentileElements[i] = (int *)calloc(nElementsInjury, sizeof(int));
    percentileTime[i] = 0.0;
    percentileValue[i] = 0.0;
  }
  PS_Old = (double *)calloc(nElementsInjury, sizeof(double));
  PSxSRArray = (double *)calloc(nElementsInjury, sizeof(double));
  elementMPS = (double *)calloc(nelements, sizeof(double));
  initialVolume = (double *)calloc(nelements, sizeof(double));

  // Array to output to Paraview
  // All threshold elements CSDM : 3, 5, 10, 15, 30, MPSR-120, MPSxSR-28
  for (int i = 0; i < threshQuantities; ++i) {
    outputDataArray[i] = thresholdElements[i];
  }
  // for (int i = 0; i < percentileQuantities; ++i) {
  //   outputDataArray[threshQuantities + i] =
  //       percentileElements[i];
  // }
  outputDataArray[threshQuantities] =
    percentileElements[0]; // Write MPS-95 to paraview
  // Percentile Values
  // Write MPS-95-Value
  outputDoubleArray[0] = PS_Old;

  // Compute the initial volume of elements for MPS file
  double volumePart[nPIDglobal];
  for (int i = 0; i < nPIDglobal; ++i) {
    volumePart[i] = 0.0;
  }
  computePartVolume(volumePart, initialVolume);
}

void CalculateInjuryCriteria(void) {
  double currentStrainMaxElem, currentStrainMinElem, currentShearMaxElem;
  double PSR = 0.0, PSxSR = 0.0;

  for (int j = 0; j < nElementsInjury; j++) {
    int i = elementIDInjury[j];
    CalculateMaximumPrincipalStrain(
        i, &currentStrainMaxElem, &currentStrainMinElem, &currentShearMaxElem);
    if (maxValue[0] < currentStrainMaxElem) {
      maxValue[0] = currentStrainMaxElem;
      maxElem[0] = i;
      maxTime[0] = Time;
    }
    if (maxValue[3] > currentStrainMinElem) {
      maxValue[3] = currentStrainMinElem;
      maxElem[3] = i;
      maxTime[3] = Time;
    }
    if (maxValue[1] < currentShearMaxElem) {
      maxValue[1] = currentShearMaxElem;
      maxElem[1] = i;
      maxTime[1] = Time;
    }
    // Update csdm element list
    for (int k = 0; k < csdmCount; ++k) {
      if (!thresholdElements[k][j]) {
        if (currentStrainMaxElem > csdmLimits[k]) {
          thresholdElements[k][j] = 1;
        }
      }
    }
    // Compute maxPSxSR
    // Compute principal strain rate using first order backaward distance
    PSR = (currentStrainMaxElem - PS_Old[j]) / dt;
    // TODO(Anil) : Sould absolute value be used for PSR ?
    PSxSR = currentStrainMaxElem * PSR;
    if (maxValue[2] < PSxSR) {
      maxValue[2] = PSxSR;
      maxElem[2] = i;
      maxTime[2] = Time;
    }
    // Calculate MPSR 120
    if (!thresholdElements[csdmCount][j]) {
      if (PSR > 120.0) {
        thresholdElements[csdmCount][j] = 1;
      }
    }
    // Calculate MPSxSR 28
    if (!thresholdElements[csdmCount + 1][j]) {
      if (PSxSR > 28.0) {
        thresholdElements[csdmCount + 1][j] = 1;
      }
    }
    PS_Old[j] = currentStrainMaxElem;
    PSxSRArray[j] = PSxSR;
    if (currentStrainMaxElem > elementMPS[i]) {
      elementMPS[i] = currentStrainMaxElem;
    }
  } // For loop over elements included for injury

  // Compute 95 percentile MPS and corresponding element list
  double MPS95 = compute95thPercentileValue(PS_Old, nElementsInjury);
  if (MPS95 > percentileValue[0]) {
    memset(percentileElements[0], 0, nElementsInjury * sizeof(int));
    percentileValue[0] = MPS95;
    percentileTime[0] = Time;
    for (int j = 0; j < nElementsInjury; j++) {
      if (PS_Old[j] >= MPS95) {
        percentileElements[0][j] = 1;
      }
    }
  }

  // Compute 95 percentile MPSxSR and corresponding element list
  double MPSxSR95 = compute95thPercentileValue(PSxSRArray, nElementsInjury);
  if (MPSxSR95 > percentileValue[1]) {
    memset(percentileElements[1], 0, nElementsInjury * sizeof(int));
    percentileValue[1] = MPSxSR95;
    percentileTime[1] = Time;
    int count = 0;
    for (int j = 0; j < nElementsInjury; j++) {
      if (PSxSRArray[j] >= MPSxSR95) {
        percentileElements[1][j] = 1;
      }
    }
  }
}

void TransformMesh(const Json::Value &jsonInput) {
  double cordAngularSensor[3];
  // Check if mesh-transformation key is present in the input json file
  // If present transform mesh, otherwise subtract center of mass
  if (!jsonInput["head-cg"].empty()) {
    cm[0] = jsonInput["head-cg"][0].asDouble();
    cm[1] = jsonInput["head-cg"][1].asDouble();
    cm[2] = jsonInput["head-cg"][2].asDouble();
    FILE_LOG_MASTER(INFO, "Center of Mass used : (%15.9e, %15.9e, %15.9e)",
                    cm[0], cm[1], cm[2]);
  } else {
    GetBodyCenterofMass(cm);
    FILE_LOG_MASTER(
        INFO, "Brain+Skull Center of Mass used : (%15.9e, %15.9e, %15.9e)",
        cm[0], cm[1], cm[2]);
  }
  if (!jsonInput["angular-sensor-position"].empty()) {
    cordAngularSensor[0] = jsonInput["angular-sensor-position"][0].asDouble();
    cordAngularSensor[1] = jsonInput["angular-sensor-position"][1].asDouble();
    cordAngularSensor[2] = jsonInput["angular-sensor-position"][2].asDouble();
    FILE_LOG_MASTER(
        INFO, "Angular sensor data about : (%15.9e, %15.9e, %15.9e)",
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
        switch (tr.at(1)) {
        case 'x':
          index[i] = 0;
          break;
        case 'y':
          index[i] = 1;
          break;
        case 'z':
          index[i] = 2;
          break;
        default:
          FILE_LOG_MASTER(ERROR, "Mesh transformation has to be x/y/z");
          TerminateFemTech(12);
          break;
        }
        switch (tr.at(0)) {
        case '-':
          factor[i] = -1;
          break;
        case '+':
          factor[i] = 1;
          break;
        default:
          FILE_LOG_MASTER(
              ERROR, "Prefix of axis not -/+, check input mesh-transformation");
          TerminateFemTech(12);
          break;
        }
      } else {
        if (tr.length() == 1) {
          factor[i] = 1;
          switch (tr.at(0)) {
          case 'x':
            index[i] = 0;
            break;
          case 'y':
            index[i] = 1;
            break;
          case 'z':
            index[i] = 2;
            break;
          default:
            FILE_LOG_MASTER(ERROR, "Mesh transformation has to be x/y/z");
            TerminateFemTech(12);
            break;
          }
        } else {
          FILE_LOG_MASTER(
              ERROR, "Error in mesh transformation. Please check input file");
          TerminateFemTech(12);
        }
      }
    }
    // Validate the input
    int indexSum = index[0] + index[1] + index[2];
    if (indexSum != 3) {
      FILE_LOG_MASTER(ERROR,
                      "Error in mesh transformation. x, y and z not present");
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
      coordTemp[transformedIndex] = cm[i] * factor[i];
    }
    // Assign cm array with transformed co-ordinates
    for (int i = 0; i < ndim; ++i) {
      cm[i] = coordTemp[i];
    }
    // Transform angular sensor position
    for (int i = 0; i < ndim; ++i) {
      int transformedIndex = index[i];
      coordTemp[transformedIndex] = cordAngularSensor[i] * factor[i];
    }
    for (int i = 0; i < ndim; ++i) {
      cordAngularSensor[i] = coordTemp[i];
    }
    // Transform all corodinates and subtract angular sensor position.
    // Transform all local co-ordinates
    for (int j = 0; j < nNodes; ++j) {
      for (int i = 0; i < ndim; ++i) {
        int transformedIndex = index[i];
        coordTemp[transformedIndex] = coordinates[i + ndim * j] * factor[i];
      }
      for (int i = 0; i < ndim; ++i) {
        coordinates[i + ndim * j] = coordTemp[i] - cordAngularSensor[i];
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
        coordinates[i + ndim * j] =
            coordinates[i + ndim * j] - cordAngularSensor[i];
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

void WriteMPS(void) {
  MPI_File mpsFilePtr;
  MPI_Info infoin;
  MPI_Info_create(&infoin);
  MPI_Info_set(infoin, "access_style", "write_once,random");
  char fileName[] = "MPSfile.dat";

  int err;
  err = MPI_File_open(MPI_COMM_WORLD, fileName, MPI_MODE_EXCL|MPI_MODE_WRONLY|MPI_MODE_CREATE, infoin, &mpsFilePtr);
  if (err != MPI_SUCCESS) {
    if (world_rank == 0) {
      MPI_File_delete(fileName, MPI_INFO_NULL);
    }
    err = MPI_File_open(MPI_COMM_WORLD, fileName, MPI_MODE_EXCL|MPI_MODE_WRONLY|MPI_MODE_CREATE, infoin, &mpsFilePtr);
    if (err != MPI_SUCCESS) {
      FILE_LOG_SINGLE(ERROR, "Unable to open MPS file");
      return;
    }
  }

  const unsigned int lineSize = 50;
  char s2[lineSize+1], s3[lineSize+1];
  unsigned int offset = 0;

  for (unsigned int i = 0; i < nelements; ++i) {
    offset = (global_eid[i]-1)*lineSize;
    sprintf(s2, "%06d, %14.5e, %14.5e\n", global_eid[i], elementMPS[i], initialVolume[i]);
    sprintf(s3, "%50s", s2);
    MPI_File_write_at(mpsFilePtr, offset, s3, lineSize, MPI_CHAR, MPI_STATUS_IGNORE);
  }

  MPI_File_close(&mpsFilePtr);
  FILE_LOG_MASTER(INFO, "MPS values written to file");
  return;
}
