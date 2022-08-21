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

using namespace std;

/*Declare Functions*/
void InitCustomPlot(const Json::Value &jsonInput);
double InitBoundaryCondition(const Json::Value &jsonInput);
void InitInjuryCriteria(void);
void CustomPlot();
void ApplyPressureBoundaryConditions();
void CalculateInjuryCriteria(void);
void WriteOutputFile();
void WriteMPS(void);
void WriteMPr(void);

/* Global Variables/Parameters */
double Time = 0.0, dt, tInitial;
int nSteps;
bool ImplicitStatic = false;
bool ImplicitDynamic = false;
bool ExplicitDynamic = true;
bool reducedIntegration = true;

double dynamicDamping = 0.000;
double ExplicitTimeStepReduction = 0.8;
double FailureTimeStep = 1e-8; // Set for max runtime of around 5 hrs on aws
double MaxTimeStep = 1e-5;
int nPlotSteps = 50;
int nWriteSteps = 100;

/* Hard coded BC nodes */
const unsigned int pressureBCNodeCount = 72;
const unsigned int fixedBCNodeCount = 117;
const unsigned int pressureBCNode[pressureBCNodeCount] = {
  5577, 2059, 5571, 2066, 2065, 2058, 5005, 5011,
  2259, 1470, 1461, 1885, 5000, 5003, 5009, 5019,
  5022, 2260, 1482, 1480, 1469, 1460, 1456, 1458,
  1884, 5002, 5001, 5004, 5010, 5018, 5020, 1479,
  1478, 1468, 1459, 1457, 4992, 5181, 5021, 5023,
  5028, 1489, 1483, 1481, 1646, 1443, 4991, 4993,
  1770, 1445, 1444, 5030, 5029, 5031, 2563, 1492,
  1491, 1490, 4989, 4980, 4979, 4988, 4981, 4990,
  1649, 1769, 1433, 1442, 1432, 1440, 1431, 1441
};
const unsigned int fixedBCNode[fixedBCNodeCount] = {
  3134, 3137, 3142, 2740, 2727, 2975, 2976, 2978,
  2980, 2983, 2985, 6414, 6413, 6412, 6411, 6410,
  6409, 6196, 6203, 6557, 6560, 6561, 6821, 6550,
  6546, 6547, 6211, 6204, 6415, 6416, 6417, 6418,
  6419, 6420, 6421, 6422, 6428, 6429, 6430, 6423,
  6424, 6431, 6432, 6433, 6425, 6426, 6434, 6427,
  6856, 6864, 6865, 6857, 6866, 6858, 3141, 3138,
  3133, 2741, 2730, 2974, 2977, 2979, 2981, 2982,
  2984, 2991, 2992, 2999, 3000, 2993, 2994, 2989,
  2990, 2988, 2995, 3001, 3002, 2996, 2987, 2986,
  2744, 2749, 3123, 3124, 3120, 2997, 3004, 3003,
  3005, 2998, 3459, 3467, 3468, 3460, 3420, 3461,
  3469, 6435, 6436, 6437, 6438, 6439, 6440, 6441,
  6869, 6870, 6871, 3474, 3473, 3472, 3012, 3011,
  3010, 3009, 3008, 3007, 3006
};

/* Global variables used only in this file */
int nodeIDtoPlot;
bool rankForCustomPlot;
bool outputRelativeDisplacement = true;
/* Global variables for bc */
int *staticBoundaryID = NULL;
int *pressureBoundaryID = NULL;
double *pressureNodalArea = NULL;
int staticBoundarySize;
int pressureBoundarySize;
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
const int maxQuantities = 5;
const std::string maxOutput[maxQuantities] = {
    "principal-max-strain", "maximum-shear-strain", "maximum-PSxSR",
    "principal-min-strain", "principal-max-pressure"};
const MPI_Op outputOperator[maxQuantities] = {MPI_MAXLOC, MPI_MAXLOC,
                                              MPI_MAXLOC, MPI_MINLOC,
                                              MPI_MAXLOC};

double maxValue[maxQuantities], maxTime[maxQuantities];
int maxElem[maxQuantities];

/* Percentile variables */
/* 0. MPS - 95 Percentile
 * 1. MPSxSR - 95 Percentile
 */
const int percentileQuantities = 3;
int *percentileElements[percentileQuantities];
double percentileTime[percentileQuantities];
double percentileValue[percentileQuantities];
const std::string percentileTag[percentileQuantities] = {"MPS-95", "MPSxSR-95", "MPr-95"};
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
double *pressureElem = NULL;
/* Store MPS, volume of each element */
double *elementMPS = NULL;
double *initialVolume = NULL;
/* Store MPr */
double *elementMPr = NULL;

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

/* Variables used to store pressure BC values */
int timeTraceBCSize;
double *timeTraceBC, *pressureTraceBC;
std::string outputFileName;

/* Variable to store MPS-95 time trace */
double *mps95TimeTrace, *pressureMPS95TimeTrace;
double MPS95;
int mpsArraySize;

int main(int argc, char **argv) {
  // Initialize FemTech including logfile and MPI
  Json::Value inputJson = InitFemTech(argc, argv);

  // Processing various solver and output flags from the json file
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
  if (!simulationJson["reduced-integration"].empty()) {
    reducedIntegration = simulationJson["reduced-integration"].asBool();
  }
  if (reducedIntegration) {
    FILE_LOG_MASTER(INFO, "Solver using reduced integration");
  }
  if (!simulationJson["dynamic-damping"].empty()) {
    dynamicDamping = simulationJson["dynamic-damping"].asDouble();
  }
  FILE_LOG_MASTER(INFO, "Dynamic damping set to : %.3f", dynamicDamping);
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
  InitCustomPlot(simulationJson);
  if (computeInjuryFlag) {
    InitInjuryCriteria();
  }

  // Initial settings for BC evaluations
  // Used if initial velocity and acceleration BC is to be set.
  tInitial = InitBoundaryCondition(simulationJson);
  Time = tInitial;

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
  if (computeInjuryFlag && (world_rank == 0)) {
    mps95TimeTrace[plot_counter] = 0.0;
    pressureMPS95TimeTrace[plot_counter] = pressureTraceBC[0];
  }

  int time_step_counter = 0;
  /** Central Difference Method - Beta and Gamma */
  // double beta = 0;
  // double gamma = 0.5;

  ShapeFunctions();
  /*  Step-1: Calculate the mass matrix similar to that of belytschko. */
  AssembleLumpedMass();
  // Needs to be after shapefunctions
  CustomPlot();

  double stableDt = StableTimeStep();
  if (stableDt < FailureTimeStep) {
    TerminateFemTech(19);
  }

  /* Obtain dt, according to Belytschko dt is calculated at end of getForce */
  dt = ExplicitTimeStepReduction * stableDt;
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
    ApplyPressureBoundaryConditions();
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
        if (computeInjuryFlag && (world_rank == 0)) {
          mps95TimeTrace[plot_counter] = MPS95;
          pressureMPS95TimeTrace[plot_counter] = interpolateLinear(\
            timeTraceBCSize, timeTraceBC, pressureTraceBC, Time);
        }
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
    stableDt = StableTimeStep();
    // Moved stable time step termination to example to write paraview output
    // file for debug
    if (stableDt < FailureTimeStep) {
      FILE_LOG_MASTER(INFO, "Time at failure : %15.6e, dt=%15.6e, tmax : %15.6e", Time, dt,
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
      TerminateFemTech(19);
    }

    dt = ExplicitTimeStepReduction * stableDt;
    time_step_counter = time_step_counter + 1;

    // Barrier not a must
    MPI_Barrier(MPI_COMM_WORLD);
  } // end explcit while loop
  // Write output if last step results not written
  int writeFileFlag = (time_step_counter - 1) % nsteps_write;
  if (writeFileFlag != 0) {
    CustomPlot();
  }
  FILE_LOG_MASTER(INFO, "End of Iterative Loop");
  mpsArraySize = plot_counter + 1;

  WriteOutputFile();
  if (computeInjuryFlag) {
    WriteMPS();
    WriteMPr();
  }

  // Free local boundary condition related arrays
  free1DArray(staticBoundaryID);
  free1DArray(timeTraceBC);
  free1DArray(pressureTraceBC);
  free1DArray(outputNodeList);
  free1DArray(outputElemList);
  // Free variables used for injury criteria
  free1DArray(PS_Old);
  free1DArray(pressureElem);
  free1DArray(elementIDInjury);
  free1DArray(PSxSRArray);
  for (int i = 0; i < threshQuantities; ++i) {
    free1DArray(thresholdElements[i]);
  }
  for (int i = 0; i < percentileQuantities; ++i) {
    free1DArray(percentileElements[i]);
  }
  free1DArray(elementMPS);
  free1DArray(elementMPr);
  free1DArray(initialVolume);
  free1DArray(pressureNodalArea);
  free1DArray(mps95TimeTrace);
  free1DArray(pressureMPS95TimeTrace);

  if (rankForCustomPlot) {
    MPI_File_close(&outputFilePtr);
  }
  FinalizeFemTech();
  return 0;
}

void ApplyPressureBoundaryConditions() {
  // Restrict stationary nodes
  for (int i = 0; i < staticBoundarySize; i++) {
    int index = staticBoundaryID[i] * ndim;
    for (int j = 0; j < ndim; ++j) {
      velocities_half[index + j] = 0.0;
    }
  }
  // Normal assumed to be in +z direction
  // Apply pressure BC using f external
  double tLocal = Time + dt;
  double localPressure = interpolateLinear(timeTraceBCSize, timeTraceBC, \
      pressureTraceBC, tLocal);
  for (int i = 0; i < pressureBoundarySize; i++) {
    int index = pressureBoundaryID[i] * ndim;
    fe[index + 2] = localPressure * pressureNodalArea[i];
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
      xCordLocal = displacements[plotNode];
      yCordLocal = displacements[plotNode + 1];
      zCordLocal = displacements[plotNode + 2];
      sprintf(output, "%s %15.9e %15.9e %15.9e", output, xCordLocal, yCordLocal, zCordLocal);
    }
    double currentStrainMaxElem, currentStrainMinElem, currentShearMaxElem;
    for (int i = 0; i < outputElemCount; ++i) {
      unsigned int plotElem = outputElemList[i];
      // Recalculation to decouple from Injury criteria
      CalculateMaximumPrincipalStrain(plotElem, &currentStrainMaxElem,
                                      &currentStrainMinElem,
                                      &currentShearMaxElem);
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
  double startTime = 0.0;
  // Read input JSON for pressure data
  if (!jsonInput["pressure"].empty()) {
    const Json::Value jsonPressure = jsonInput["pressure"];
    if (!jsonPressure["time"].empty()) {
      timeTraceBCSize = jsonPressure["time"].size();
    } else{
      FILE_LOG_SINGLE(ERROR, "Pressure-time array absent");
      TerminateFemTech(12);
    }
    timeTraceBC = (double*)malloc(sizeof(double) * timeTraceBCSize);
    jsonToArray(timeTraceBC, jsonPressure["time"]);

    pressureTraceBC = (double*)malloc(sizeof(double) * timeTraceBCSize);
    if (!jsonPressure["average"].empty()) {
      unsigned int pressureTraceSize = 0;
      pressureTraceSize = jsonPressure["average"].size();
      assert(pressureTraceSize == timeTraceBCSize);

      jsonToArray(pressureTraceBC, jsonPressure["average"]);
    } else {
      if (!jsonPressure["head"].empty()) {
        unsigned int pressureTrace1Size = 0;
        pressureTrace1Size = jsonPressure["head"].size();
        assert(pressureTrace1Size == timeTraceBCSize);
      } else{
        FILE_LOG_SINGLE(ERROR, "Head pressure trace absent");
        TerminateFemTech(12);
      }
      if (!jsonPressure["shoulder"].empty()) {
        unsigned int pressureTrace2Size = 0;
        pressureTrace2Size = jsonPressure["shoulder"].size();
        assert(pressureTrace2Size == timeTraceBCSize);
      } else{
        FILE_LOG_SINGLE(ERROR, "Shoulder pressure trace absent");
        TerminateFemTech(12);
      }
      if (!jsonPressure["chest"].empty()) {
        unsigned int pressureTrace3Size = 0;
        pressureTrace3Size = jsonPressure["chest"].size();
        assert(pressureTrace3Size == timeTraceBCSize);
      } else{
        FILE_LOG_SINGLE(ERROR, "Chest pressure trace absent");
        TerminateFemTech(12);
      }
      double *pressureTrace1 = (double*)malloc(sizeof(double) * timeTraceBCSize);
      double *pressureTrace2 = (double*)malloc(sizeof(double) * timeTraceBCSize);
      double *pressureTrace3 = (double*)malloc(sizeof(double) * timeTraceBCSize);

      jsonToArray(pressureTrace1, jsonPressure["head"]);
      jsonToArray(pressureTrace2, jsonPressure["shoulder"]);
      jsonToArray(pressureTrace3, jsonPressure["chest"]);

      // Average over the three pressure sensor values to arrive at pressure to be
      // applied on the fore-head
      for (unsigned int i = 0; i < timeTraceBCSize; ++i) {
        pressureTraceBC[i] = (pressureTrace1[i] + pressureTrace2[i] + pressureTrace3[i])/3.0;
      }

      free(pressureTrace1);
      free(pressureTrace2);
      free(pressureTrace3);
    }

    // Convert time from milli-seconds to seconds
    for (int i = 0; i < timeTraceBCSize; ++i) {
      timeTraceBC[i] = 0.001 * timeTraceBC[i];
    }
    startTime = timeTraceBC[0];
  } else { 
    FILE_LOG_MASTER(ERROR, "Pressure BC time traces absent in the input file");
    TerminateFemTech(12);
  }

  // Find count of nodes with specified partID on the current process
  staticBoundarySize = 0;
  // One time operation, can be optimized if both arrays are sorted
  for (unsigned int i = 0; i < fixedBCNodeCount; ++i) {
    const unsigned int nodeItem = fixedBCNode[i];
    for (int j = 0; j < nNodes; ++j) {
      if (globalNodeID[j] == nodeItem) {
        ++staticBoundarySize;
      }
    }
  }
  // Allocate memory to store node IDs of staionary rigid nodes.
  // Loop over this array in apply pressure BC function
  staticBoundaryID = (int *)malloc(staticBoundarySize * sizeof(int));
  if (staticBoundaryID == NULL) {
    FILE_LOG_SINGLE(ERROR, "Unable to alocate staticBoundaryID");
    TerminateFemTech(12);
  }
  // rigidnodecount will be <= fixedBCNodeCount due to the parallel mesh
  // partitioning
  // Store all nodes to be made rigid
  // set boundary array value to 1, so that displacements are not updated for
  // displacement prescribed nodes
  int nodePtr = 0;
  for (unsigned int i = 0; i < fixedBCNodeCount; ++i) {
    const unsigned int nodeItem = fixedBCNode[i];
    for (int j = 0; j < nNodes; ++j) {
      if (globalNodeID[j] == nodeItem) {
        staticBoundaryID[nodePtr] = j; // Store all rigid nodes
        nodePtr = nodePtr + 1;
        int index = j * ndim;
        boundary[index] = 1;
        boundary[index + 1] = 1;
        boundary[index + 2] = 1;
      }
    }
  }
  assert(nodePtr == staticBoundarySize);

  FILE_LOG(INFO, "%d nodes given rigid motion", staticBoundarySize);

  // Set pressure BC nodes
  // Find count of nodes with specified pressure on the current process
  pressureBoundarySize = 0;
  // One time operation, can be optimized if both arrays are sorted
  for (unsigned int i = 0; i < pressureBCNodeCount; ++i) {
    const unsigned int nodeItem = pressureBCNode[i];
    for (int j = 0; j < nNodes; ++j) {
      if (globalNodeID[j] == nodeItem) {
        ++pressureBoundarySize;
      }
    }
  }
  // Allocate memory to store node IDs of nodes with prescibed pressure.
  // Loop over this array in apply pressure BC function
  pressureBoundaryID = (int *)malloc(pressureBoundarySize * sizeof(int));
  if (pressureBoundaryID == NULL) {
    FILE_LOG_SINGLE(ERROR, "Unable to alocate pressureBoundaryID");
    TerminateFemTech(12);
  }
  // pressure boundary size will be <= pressureBCNodeCount due to the parallel mesh
  // partitioning
  // Store all nodes to prescribe pressure BC
  nodePtr = 0;
  for (unsigned int i = 0; i < pressureBCNodeCount; ++i) {
    const unsigned int nodeItem = pressureBCNode[i];
    for (int j = 0; j < nNodes; ++j) {
      if (globalNodeID[j] == nodeItem) {
        pressureBoundaryID[nodePtr] = j; // Store all rigid nodes
        nodePtr = nodePtr + 1;
      }
    }
  }
  assert(nodePtr == pressureBoundarySize);

  FILE_LOG(INFO, "%d nodes given pressure BC", pressureBoundarySize);
  // Split the MPI communicator into those containing pressure boundary nodes
  // int splitColor = (int)(pressureBoundarySize > 0);
  // MPI_Comm pressureBC_Comm;
  // MPI_Comm_split(MPI_COMM_WORLD, splitColor, world_rank, &pressureBC_Comm);
  // // Work on procs with pressure BC nodes only
  // if (splitColor) {
  // }

  int pressureNodeToSend = pressureBoundarySize;
  int *pressureBCSize = NULL, *cumPressureBCSize = NULL, *pressureGlobalID = NULL;
  int *pressureBCCoordSize = NULL, *cumPressureBCCoordSize = NULL;
  double *pressureNodeCoord = NULL, *areaOnMaster = NULL;
  if (world_rank == 0) {
    pressureBCSize = (int*)malloc(sizeof(int)*world_size);
    cumPressureBCSize = (int*)malloc(sizeof(int)*(world_size+1));
    pressureBCCoordSize = (int*)malloc(sizeof(int)*world_size);
    cumPressureBCCoordSize = (int*)malloc(sizeof(int)*(world_size+1));
  }
  // Get size of pressure BC points on each process on master
  MPI_Gather(&pressureBoundarySize, 1, MPI_INT, pressureBCSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
  // Create a MPI Communicator with only processes with pressure boundary nodes

  int totalPressureBCPoints = 0;
  if (world_rank == 0) {
    // Compute the cumulative number of points to act as pointers to the array 
    cumPressureBCSize[0] = 0;
    cumPressureBCCoordSize[0] = 0;
    for (unsigned int i = 0; i < world_size; ++i) {
      cumPressureBCSize[i+1] = cumPressureBCSize[i] + pressureBCSize[i];
      cumPressureBCCoordSize[i+1] = ndim * cumPressureBCSize[i+1];
      pressureBCCoordSize[i] = ndim * pressureBCSize[i];
    }
    // Use the last element of cumulative array to find total count of nodes.
    // This can be greater than number of nodes on which pressure is applied as
    // duplicate nodes can be present.
    totalPressureBCPoints = cumPressureBCSize[world_size];
    pressureGlobalID = (int*)malloc(sizeof(int)*totalPressureBCPoints);
    pressureNodeCoord = (double*)malloc(sizeof(double)*totalPressureBCPoints*ndim);
    areaOnMaster = (double*)malloc(sizeof(double)*totalPressureBCPoints);
  }
  int *globalIDtoMaster;
  double *coordToMaster;
  globalIDtoMaster = (int*)malloc(sizeof(int)*pressureBoundarySize);
  coordToMaster = (double*)malloc(sizeof(double)*pressureBoundarySize*ndim);
  // Copy global element ID as well as it co-ordinates to arrays to send to
  // master node to compute area
  for (unsigned int i = 0; i < pressureBoundarySize; ++i) {
    const unsigned int nodeID = pressureBoundaryID[i];
    globalIDtoMaster[i] = globalNodeID[nodeID];
    const unsigned int index0 = nodeID * ndim;
    const unsigned int index1 = i * ndim;
    for (unsigned int j = 0; j < ndim; ++j) {
      coordToMaster[index1 + j] = coordinates[index0 + j];
    }
  }
  // FILE_LOGArrayIntSingle(WARNING, globalIDtoMaster, pressureBoundarySize, "ID to master : ");
  // Get global nodeIDs to master from all processes
  MPI_Gatherv(globalIDtoMaster, pressureBoundarySize, MPI_INT, \
      pressureGlobalID, pressureBCSize, cumPressureBCSize, MPI_INT, \
      0, MPI_COMM_WORLD);
  // FILE_LOGArrayIntMaster(WARNING, pressureGlobalID, totalPressureBCPoints, "Array on master : ");
  // Get global corresponding co-ordinates to master from all processes
  MPI_Gatherv(coordToMaster, ndim*pressureBoundarySize, MPI_DOUBLE, \
      pressureNodeCoord, pressureBCCoordSize, cumPressureBCCoordSize, \
      MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (world_rank == 0) {
    const double piBy4 = atan(1.0);
    // Find the nearest point to a given point in the list using a O(N^2) 
    // brute force algorithm. We are dealing with small arrays here and its
    // a one time operation.
    for (unsigned int i = 0; i < totalPressureBCPoints; ++i) {
      double minDist = 1e12;
      const unsigned int index0 = i * ndim;
      for (unsigned int j = 0; j < totalPressureBCPoints; ++j) {
        if (i == j) {
          continue;
        }
        if (pressureGlobalID[i] == pressureGlobalID[j]) {
          continue;
        }
        double distance = 0.0;
        const unsigned int index1 = j * ndim;
        for (unsigned int k = 0; k < ndim; ++k) {
          distance += pow(pressureNodeCoord[index0 + k] - \
              pressureNodeCoord[index1 + k], 2);
        }
        if (distance < minDist) {
          minDist = distance;
        }
      }
      areaOnMaster[i] = piBy4*minDist;
    }
  }
  // FILE_LOGArrayMaster(WARNING, areaOnMaster, totalPressureBCPoints, "Area on master : ");
  pressureNodalArea = (double*)malloc(sizeof(double)*pressureBoundarySize);
  // Scatter all nodal areas to respective processes 
  MPI_Scatterv(areaOnMaster, pressureBCSize, cumPressureBCSize, MPI_DOUBLE, \
      pressureNodalArea, pressureBoundarySize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  free1DArray(globalIDtoMaster);
  free1DArray(coordToMaster);
  if (world_rank == 0) {
    free1DArray(pressureBCSize);
    free1DArray(cumPressureBCSize);
    free1DArray(pressureGlobalID);
    free1DArray(pressureNodeCoord);
    free1DArray(pressureBCCoordSize);
    free1DArray(cumPressureBCCoordSize);
    free1DArray(areaOnMaster);
  }
  return startTime;
}

void WriteOutputFile() {
  Json::Value output;
  Json::Value vec(Json::arrayValue);

  // Write basic info that does not require any computation
  if (world_rank == 0) {
    output["output-file"] = outputFileName + ".pvd";

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
    double totalVolume = 0.0;
    if (world_rank == 0) {
      // Excluded PID hardcoded for CSDM computation
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
    const int outputJsonCount = csdmCount;
    std::string jsonOutputTag[outputJsonCount];
    int *outputJsonArray[outputJsonCount];
    for (int i = 0; i < csdmCount; ++i) {
      jsonOutputTag[i] = thresholdTag[i];
      outputJsonArray[i] = thresholdElements[i];
    }

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
        Json::Value outputThreshold;
        outputThreshold["global-element-id"] = elementList;
        Json::StreamWriterBuilder builder;
        builder["commentStyle"] = "None";
        builder["indentation"] = "  ";
        builder.settings_["precision"] = 4;
        std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
        std::ofstream oFile(jsonOutputTag[i] + ".json");

        if (i < csdmCount) {
          outputThreshold["value"] = thresholdVolume[i] / totalVolume;
        } else {
          int ind = i - csdmCount;
          outputThreshold["value"] = percentileValue[ind];
          outputThreshold["time"] = percentileTime[ind];
        }
        writer->write(outputThreshold, &oFile);
      }
      free(localElemList);
      MPI_Barrier(MPI_COMM_WORLD);
    }
    if (world_rank == 0) {
      free(countPerProc);
      free(cumCountPerProc);
    }
    // Write MPS95 value instead of MPS value for better post processing
    output[maxOutput[0]]["value"] = percentileValue[0];
  }

  // Write output to file
  if (world_rank == 0) {
    Json::StreamWriterBuilder builder;
    builder["commentStyle"] = "None";
    builder["indentation"] = "  ";
    builder.settings_["precision"] = 4;
    std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
    std::ofstream oFile(uid + "_output.json");
    writer->write(output, &oFile);
  }
  FILE_LOG_MASTER(INFO, "Json output files written");

  // Write basic info that does not require any computation
  if (world_rank == 0) {
    Json::Value mps95TimeTraceJson;
    Json::Value vecTime(Json::arrayValue);
    Json::Value vecMPS95(Json::arrayValue);
    Json::Value vecPressure(Json::arrayValue);

    vecTime.resize(mpsArraySize);
    vecMPS95.resize(mpsArraySize);
    vecPressure.resize(mpsArraySize);
    for (int i = 0; i < mpsArraySize; ++i) {
      vecTime[i] = stepTime[i];
      vecMPS95[i] = mps95TimeTrace[i];
      vecPressure[i] = pressureMPS95TimeTrace[i];
    }

    mps95TimeTraceJson["time"] = vecTime;
    mps95TimeTraceJson["MPS95"] = vecMPS95;
    mps95TimeTraceJson["pressure"] = vecPressure;

    Json::StreamWriterBuilder builder;
    builder["commentStyle"] = "None";
    builder["indentation"] = "  ";
    builder.settings_["precision"] = 4;
    std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
    std::ofstream oFile(uid + "_mpsTimeTrace.json");
    writer->write(mps95TimeTraceJson, &oFile);
  }
  FILE_LOG_MASTER(INFO, "MPS time trace written as Json file");
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
    // Using element list only for MPS-95, to write to paraview
    if (i == 0) {
      percentileElements[i] = (int *)calloc(nElementsInjury, sizeof(int));
    }
    percentileTime[i] = 0.0;
    percentileValue[i] = 0.0;
  }
  PS_Old = (double *)calloc(nElementsInjury, sizeof(double));
  PSxSRArray = (double *)calloc(nElementsInjury, sizeof(double));
  pressureElem = (double *)calloc(nElementsInjury, sizeof(double));
  elementMPS = (double *)calloc(nelements, sizeof(double));
  elementMPr = (double *)calloc(nelements, sizeof(double));
  initialVolume = (double *)calloc(nelements, sizeof(double));

  // Array to output to Paraview
  // All threshold elements CSDM : 3, 5, 10, 15, 30, MPSR-120, MPSxSR-28
  for (int i = 0; i < threshQuantities; ++i) {
    outputDataArray[i] = thresholdElements[i];
  }
  outputDataArray[threshQuantities] = percentileElements[0]; // Write MPS-95 to paraview
  // Percentile Values
  // Write MPS-95-Value
  outputDoubleArray[0] = PS_Old;

  // Compute the initial volume of elements for MPS file
  double volumePart[nPIDglobal];
  for (int i = 0; i < nPIDglobal; ++i) {
    volumePart[i] = 0.0;
  }
  computePartVolume(volumePart, initialVolume);
  // Allocate variables for MPS95 time trace
  if (world_rank == 0) {
    mps95TimeTrace = (double*)malloc(MAXPLOTSTEPS*sizeof(double));
    pressureMPS95TimeTrace = (double*)malloc(MAXPLOTSTEPS*sizeof(double));
  }
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
    // TODO(Anil) : Should absolute value be used for PSR ?
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
    // compute and store element pressure quantities
    double elementPressure = CalculateElementPressure(i);
    pressureElem[j] = elementPressure;
    if (elementPressure > elementMPr[i]) {
      elementMPr[i] = elementPressure;
    }
    if (maxValue[4] < elementPressure) {
      maxValue[4] = elementPressure;
      maxElem[4] = i;
      maxTime[4] = Time;
    }
  } // For loop over elements included for injury

  // Compute 95 percentile MPS and corresponding element list
  MPS95 = compute95thPercentileValue(PS_Old, nElementsInjury);
  if (MPS95 > percentileValue[0]) {
    percentileValue[0] = MPS95;
    percentileTime[0] = Time;
  }
  // Store instantaneous elements with MPS > MPS95
  // Global MPS95 element list is computed from MPS text file
  // Used to plot in paraview
  memset(percentileElements[0], 0, nElementsInjury * sizeof(int));
  for (int j = 0; j < nElementsInjury; j++) {
    if (PS_Old[j] >= MPS95) {
      percentileElements[0][j] = 1;
    }
  }
  // Compute 95 percentile maximum pressure 
  double MPr95 = compute95thPercentileValue(pressureElem, nElementsInjury);
  if (MPr95 > percentileValue[2]) {
    percentileValue[2] = MPr95;
    percentileTime[2] = Time;
  }

  // Compute 95 percentile MPSxSR and corresponding element list
  double MPSxSR95 = compute95thPercentileValue(PSxSRArray, nElementsInjury);
  if (MPSxSR95 > percentileValue[1]) {
    // memset(percentileElements[1], 0, nElementsInjury * sizeof(int));
    percentileValue[1] = MPSxSR95;
    percentileTime[1] = Time;
    // TODO(Anil) MPSxSR95 may only be representitive
    // for (int j = 0; j < nElementsInjury; j++) {
    //   if (PSxSRArray[j] >= MPSxSR95) {
    //     percentileElements[1][j] = 1;
    //   }
    // }
  }
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

void WriteMPr(void) {
  MPI_File mpsFilePtr;
  MPI_Info infoin;
  MPI_Info_create(&infoin);
  MPI_Info_set(infoin, "access_style", "write_once,random");
  char fileName[] = "MPrFile.dat";

  int err;
  err = MPI_File_open(MPI_COMM_WORLD, fileName, MPI_MODE_EXCL|MPI_MODE_WRONLY|MPI_MODE_CREATE, infoin, &mpsFilePtr);
  if (err != MPI_SUCCESS) {
    if (world_rank == 0) {
      MPI_File_delete(fileName, MPI_INFO_NULL);
    }
    err = MPI_File_open(MPI_COMM_WORLD, fileName, MPI_MODE_EXCL|MPI_MODE_WRONLY|MPI_MODE_CREATE, infoin, &mpsFilePtr);
    if (err != MPI_SUCCESS) {
      FILE_LOG_SINGLE(ERROR, "Unable to open MPr file");
      return;
    }
  }

  const unsigned int lineSize = 50;
  char s2[lineSize+1], s3[lineSize+1];
  unsigned int offset = 0;

  for (unsigned int i = 0; i < nelements; ++i) {
    offset = (global_eid[i]-1)*lineSize;
    sprintf(s2, "%06d, %14.5e, %14.5e\n", global_eid[i], elementMPr[i], initialVolume[i]);
    sprintf(s3, "%50s", s2);
    MPI_File_write_at(mpsFilePtr, offset, s3, lineSize, MPI_CHAR, MPI_STATUS_IGNORE);
  }

  MPI_File_close(&mpsFilePtr);
  FILE_LOG_MASTER(INFO, "MPr values written to file");
  return;
}
