#include "FemTech.h"

#include <assert.h>
// #include <fenv.h>

/*Delare Functions*/
void ApplyVelocityBoundaryConditions(double);
void InitVelocityBoundaryConditions();
void CustomPlot();

double Time, dt;
int nSteps;

double dynamicDamping = 0.000;
double ExplicitTimeStepReduction = 0.668;
double FailureTimeStep = 1e-11;
double MaxTimeStep = 1e-1;

int nPlotSteps = 50;
int nFieldSkip = 1; // 1 in nFieldSkip plot steps will be used to output VTU
bool ImplicitStatic = false;
bool ImplicitDynamic = false;
bool ExplicitDynamic = true;
bool RI = true;

// Parameters of simple tension test
double tMax = 1.00;  // max simulation time in seconds
double dMax = 0.001; // max displacment in meters

int main(int argc, char **argv) {
  // feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  InitFemTechWoInput(argc, argv);

  ReadInputFile(argv[1]);
  ReadMaterials();

  PartitionMesh();

  AllocateArrays();

  std::string meshFile(argv[1]);
  size_t lastindex = meshFile.find_last_of(".");
  std::string outputFileName = meshFile.substr(0, lastindex);
  /* Write inital, undeformed configuration*/
  Time = 0.0;
  int plot_counter = 0;
  WriteVTU(outputFileName.c_str(), plot_counter);
  stepTime[plot_counter] = Time;
  CustomPlot();

  // Dynamic Explcit solution Belytschko Box 6.1
  dt = 0.0;

  // Initial operations 
  int time_step_counter = 0;
  ShapeFunctions();
  /*  Step-1: Calculate the mass matrix similar to that of belytschko. */
  AssembleLumpedMass();
  InitVelocityBoundaryConditions();

  /* Obtain dt, according to Belytschko dt is calculated at end of getForce */
  dt = ExplicitTimeStepReduction * StableTimeStep();

  // Used if initial velocity BC is to be set.
  ApplyVelocityBoundaryConditions(0.5*dt);

  // For cases with no residual stress and external forces, 
  // GetFoce is not required, we can directly set accelerations at zero time
  // step to be zero
  /* Step-2: getforce step from Belytschko */
  // GetForce(); // Calculating the force term.

  /* Step-3: Calculate accelerations */
  // CalculateAccelerations();
  for (int i = 0; i < nDOF; i++) {
    accelerations[i] = 0.0;
  }

  nSteps = (int)(tMax / dt);
  double dtPlot = tMax/double(nPlotSteps);
  double nextPlotTime = dtPlot;
  int nsteps_plot = (int)(nSteps / nPlotSteps);

  double t_n = 0.0;

  FILE_LOG_MASTER(INFO, "---------------------------------");
  FILE_LOG_MASTER(INFO, "Tmax : %15.6e, Initial dt : %15.6e", tMax, dt);
  FILE_LOG_MASTER(INFO, "nSteps = %d, nsteps_plot = %d", nSteps, nsteps_plot);
  FILE_LOG_MASTER(INFO, "-------------- Loop -------------");

  bool writeFlag = false;
  int  vtuSkipCount = 0;

  /* Step-4: Time loop starts....*/
  while (Time < tMax) {
    t_n = Time;
    double t_np1 = Time + dt;
    if (t_np1 > nextPlotTime) {
      // Adjust dt to plot at specific time steps
      t_np1 = nextPlotTime;
      dt = nextPlotTime - Time;
      // Update next plot time step
      nextPlotTime = nextPlotTime + dtPlot;
      writeFlag = true;
    }
    Time = t_np1; /*Update the time by adding full time step */
    FILE_LOG_MASTER(INFO, "Time : %15.6e, dt=%15.6e, tmax : %15.6e", Time, dt,
                    tMax);
    double dtby2 = 0.5*dt;
    double t_nphalf = t_n + dtby2; // equ 6.2.1

    /* Step 6 Enforce velocity boundary Conditions */
    ApplyVelocityBoundaryConditions(t_nphalf);

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
    if(RI)
	memcpy(f_hgprev, f_hg, nDOF * sizeof(double));

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
    CheckEnergy(Time, writeFlag);
    
    // int writeFlag = time_step_counter % nsteps_plot;
    if (writeFlag) {
      CustomPlot();
      vtuSkipCount = vtuSkipCount + 1;
      if (!(vtuSkipCount%nFieldSkip)) {
        plot_counter = plot_counter + 1;
        CalculateStrain();
        FILE_LOG(INFO, "------ Plot %d: WriteVTU", plot_counter);
        WriteVTU(outputFileName.c_str(), plot_counter);
        if (plot_counter < MAXPLOTSTEPS) {
          stepTime[plot_counter] = Time;
          WritePVD(outputFileName.c_str(), plot_counter);
        }
        vtuSkipCount = 0;
      }

      FILE_LOGMatrixRM(DEBUGLOG, displacements, nNodes, ndim,
                        "Displacement Solution");
      writeFlag = false;
    }
    // time_step_counter = time_step_counter + 1;
    dt = ExplicitTimeStepReduction * StableTimeStep();
    // Barrier not a must
    MPI_Barrier(MPI_COMM_WORLD);
  } // end explcit while loop
  FILE_LOG_MASTER(INFO, "End of Iterative Loop");

  // Write out the last time step
  // CustomPlot();
  FILE_LOGMatrixRM(DEBUGLOG, displacements, nNodes, ndim,
                   "Final Displacement Solution");

  FinalizeFemTech();
  return 0;
}

void ApplyVelocityBoundaryConditions(double) {
  double tol = 1e-5;
  int index;

  for (int i = 0; i < nNodes; i++) {
    // if x value = 0, constrain node to x plane (0-direction)
    index = ndim * i + 1;
    if (fabs(coordinates[index] - 0.0) < tol) {
      velocities_half[index-1] = 0.0;
    }
    // if y coordinate = 0, constrain node to y plane (1-direction)
    index = ndim * i + 1;
    if (fabs(coordinates[index] - 0.0) < tol) {
      velocities_half[index] = 0.0;
    }
    // if z coordinate = 0, constrain node to z plane (2-direction)
    index = ndim * i + 2;
    if (fabs(coordinates[index] - 0.0) < tol) {
      velocities_half[index] = 0.0;
    }
    if (fabs(coordinates[index] - 0.005) < tol) {
      velocities_half[index] = 0.0;
    }
    // if y coordinate = 1, apply disp. to node = 0.1 (1-direction)
    // if y coordinate = 1, apply disp. to node = 0.1 (1-direction)
    index = ndim * i + 1;
    if (fabs(coordinates[index] - 0.001) < tol) {
      velocities_half[index-1] = dMax / tMax;
      velocities_half[index] = 0.0;
    }
  }
  FILE_LOG_MASTER(INFO, "Time = %10.5E, Applied Velocity = %10.5E", Time,
                  dMax / tMax);
  return;
}

void InitVelocityBoundaryConditions() {
  double tol = 1e-5;
  int index;
  for (int i = 0; i < nNodes; i++) {
    // if x value = 0, constrain node to x plane (0-direction)
    index = ndim * i + 1;
    if (fabs(coordinates[index] - 0.0) < tol) {
      boundary[index-1] = 1;
    }
    // if y coordinate = 0, constrain node to y plane (1-direction)
    index = ndim * i + 1;
    if (fabs(coordinates[index] - 0.0) < tol) {
      boundary[index] = 1;
    }
    // if z coordinate = 0, constrain node to z plane (2-direction)
    index = ndim * i + 2;
    if (fabs(coordinates[index] - 0.0) < tol) {
      boundary[index] = 1;
    }
    if (fabs(coordinates[index] - 0.005) < tol) {
      boundary[index] = 1;
    }
    // if y coordinate = 1, apply disp. to node = 0.1 (1-direction)
    index = ndim * i + 1;
    if (fabs(coordinates[index] - 0.001) < tol) {
      boundary[index-1] = 1;
      boundary[index] = 1;
    }
  }
  return;
}

void CustomPlot() {
  double tol = 1e-5;
  FILE *datFile;
  int x = 0;
  int y = 1;
  int z = 2;

  if (fabs(Time - 0.0) < 1e-16) {
    datFile = fopen("plot.dat", "w");
    fprintf(datFile, "# Results for Node ?\n");
    fprintf(datFile, "# Time  DispX  DispY  DispZ  SigmaXX  SigmaYY  SigmaZZ  "
                     "SigmaXY  Sigma XZ  SigmaYZ\n");
    fprintf(datFile, "%13.5e  %13.5e  %13.5e  %13.5e  %13.5e  %13.5e  %13.5e  "
            "%13.5e  %13.5e  %13.5e\n", 
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  } else {
    datFile = fopen("plot.dat", "a");
    for (int i = 0; i < nNodes; i++) {
      if (fabs(coordinates[ndim * i + x] - 0.002) < tol &&
          fabs(coordinates[ndim * i + y] - 0.0005) < tol &&
          fabs(coordinates[ndim * i + z] - 0.002) < tol) {

        fprintf(datFile, "%13.5e  %13.5e  %13.5e  %13.5e  ", Time,
                displacements[ndim * i + x], displacements[ndim * i + y],
                displacements[ndim * i + z]);
      }
    }
    // Compute Cauchy stress for the element
    double stressElem[6];
    CalculateElementStress(94, stressElem);
    // Print all six stress components of 1st element stress
    fprintf(datFile, "%13.5e  %13.5e  %13.5e  %13.5e  %13.5e  %13.5e\n",
            stressElem[0], stressElem[1], stressElem[2], stressElem[5],
            stressElem[4], stressElem[3]);
  }
  fclose(datFile);
  return;
}
