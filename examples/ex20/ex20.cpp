#include "FemTech.h"

#include <assert.h>

/*Delare Functions*/
void ApplyBoundaryConditions(double tMax);
void CustomPlot();

double Time, dt;
int nSteps;
double ExplicitTimeStepReduction = 0.7;
double FailureTimeStep = 1e-11;

int nPlotSteps = 50;
bool ImplicitStatic = false;
bool ImplicitDynamic = false;
bool ExplicitDynamic = true;

/*Drupal*/
bool reducedIntegration = true;
double tInitial = 0.000, dynamicDamping = 0.000;

double MaxTimeStep = 1e-5;

int main(int argc, char **argv) {
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

  // Dynamic Explcit solution using....

  dt = 0.0;
  double tMax = 1.0; // max simulation time in seconds

  int time_step_counter = 0;
  /** Central Difference Method - Beta and Gamma */
  // double beta = 0;
  // double gamma = 0.5;

  ShapeFunctions();
  /*  Step-1: Calculate the mass matrix similar to that of belytschko. */
  AssembleLumpedMass();

  // Used if initial velocity and acceleration BC is to be set.
  ApplyBoundaryConditions(tMax);

  /* Obtain dt, according to Belytschko dt is calculated at end of getForce */
  dt = ExplicitTimeStepReduction * StableTimeStep();
  /* Step-2: getforce step from Belytschko */
  GetForce(); // Calculating the force term.

  /* Step-3: Calculate accelerations */
  CalculateAccelerations();

  nSteps = (int)(tMax / dt);
  int nsteps_plot = (int)(nSteps / nPlotSteps);
  if (nsteps_plot == 0) {
    nsteps_plot = 1;
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
    FILE_LOG_MASTER(INFO, "Time : %15.6e, dt=%15.6e, tmax : %15.6e", Time, dt, tMax);
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
    /* Step 6 Enforce displacement boundary Conditions */
    ApplyBoundaryConditions(tMax);

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
      plot_counter = plot_counter + 1;
      CalculateStrain();
      FILE_LOG(INFO, "------ Plot %d: WriteVTU", plot_counter);
      WriteVTU(outputFileName.c_str(), plot_counter);
      if (plot_counter < MAXPLOTSTEPS) {
        stepTime[plot_counter] = Time;
        WritePVD(outputFileName.c_str(), plot_counter);
      }
      CustomPlot();

      FILE_LOGMatrixRM(DEBUGLOG, displacements, nNodes, ndim, "Displacement Solution");
    }
    time_step_counter = time_step_counter + 1;
    dt = ExplicitTimeStepReduction * StableTimeStep();
    // Barrier not a must
    MPI_Barrier(MPI_COMM_WORLD);
  } // end explcit while loop
  FILE_LOG_MASTER(INFO, "End of Iterative Loop");

  // Write out the last time step
  CustomPlot();
  FILE_LOGMatrixRM(DEBUGLOG, displacements, nNodes, ndim, "Final Displacement Solution");

  FinalizeFemTech();
  return 0;
}

void ApplyBoundaryConditions(double tMax) {
  double tol = 1e-5;
  double load = 625.0;
  // Apply Ramped Pressure
  double rampedUnit = Time*load/tMax;

  int index;
  for (int i = 0; i < nNodes; i++) {
    // if x value = 0, constrain node to x plane (0-direction)
    index = ndim * i + 0;
    if (fabs(coordinates[index] - 0.0) < tol) {
      boundary[index] = 1;
      displacements[index] = 0.0;
      velocities[index] = 0.0;
      // For energy computations
      accelerations[index] = 0.0;
    }
    // if y coordinate = 0, constrain node to y plane (1-direction)
    index = ndim * i + 1;
    if (fabs(coordinates[index] - 0.0) < tol) {
      boundary[index] = 1;
      displacements[index] = 0.0;
      velocities[index] = 0.0;
      accelerations[index] = 0.0;
    }
    // if z coordinate = 0, constrain node to z plane (2-direction)
    index = ndim * i + 2;
    if (fabs(coordinates[index] - 0.0) < tol) {
      boundary[index] = 1;
      displacements[index] = 0.0;
      velocities[index] = 0.0;
      accelerations[index] = 0.0;
    }
    // if y coordinate = 1, apply force in y direction
    index = ndim * i + 1;
    if (fabs(coordinates[index] - 0.005) < tol) {
      fe[index] = rampedUnit;
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
    fprintf(datFile, "# Time  DispX    DispY   DispZ\n");
    fprintf(datFile, "%11.5e %11.5e  %11.5e  %11.5e\n", 0.0, 0.0, 0.0, 0.0);

  } else {
    datFile = fopen("plot.dat", "a");
    for (int i = 0; i < nNodes; i++) {
      if (fabs(coordinates[ndim * i + x] - 0.005) < tol &&
          fabs(coordinates[ndim * i + y] - 0.005) < tol &&
          fabs(coordinates[ndim * i + z] - 0.005) < tol) {

        fprintf(datFile, "%11.5e %11.5e  %11.5e  %11.5e\n", Time,
                displacements[ndim * i + x], displacements[ndim * i + y],
                displacements[ndim * i + z]);
      }
    }
  }

  fclose(datFile);
  return;
}
