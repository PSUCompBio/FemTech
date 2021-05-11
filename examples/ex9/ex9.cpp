#include "FemTech.h"

#include <assert.h>
#include <fenv.h>

/*Delare Functions*/
void ApplyBoundaryConditions();
void ApplyVelocityBoundaryConditions(double);
void InitVelocityBoundaryConditions();
void CustomPlot();

double Time, dt;
int nSteps;

double dynamicDamping = 0.001;
double ExplicitTimeStepReduction = 0.7;
double FailureTimeStep = 1e-11;
double MaxTimeStep = 1e-1;

int nPlotSteps = 1000;
bool ImplicitStatic = false;
bool ImplicitDynamic = false;
bool ExplicitDynamic = true;

// Parameters of simple tension test
double tMax = 1.00;  // max simulation time in seconds
double dMax = 0.007; // max displacment in meters

int main(int argc, char **argv) {
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
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

  if (ImplicitStatic) {
    // Static solution
    ShapeFunctions();
    CreateLinearElasticityCMatrix();
    ApplyBoundaryConditions();
    Assembly((char *)"stiffness");
    ApplySteadyBoundaryConditions();
    SolveSteadyImplicit();
    Time = tMax;
    /* Write final, deformed configuration*/
    WriteVTU(argv[1], 1);
  } else if (ImplicitDynamic) {
    // Dynamic Implicit solution using Newmark's scheme for time integration
    dt = 0.1;
    ShapeFunctions();
    CreateLinearElasticityCMatrix();
    Time = 1.0;
    ApplyBoundaryConditions();
    Assembly((char *)"stiffness");
    Assembly((char *)"mass");
    /* beta and gamma of Newmark's scheme */
    double beta = 0.25;
    double gamma = 0.5;
    SolveUnsteadyNewmarkImplicit(beta, gamma, dt, tMax, argv[1]);
  } else if (ExplicitDynamic) {
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
        plot_counter = plot_counter + 1;
        CalculateStrain();
        FILE_LOG(INFO, "------ Plot %d: WriteVTU", plot_counter);
        WriteVTU(outputFileName.c_str(), plot_counter);
        if (plot_counter < MAXPLOTSTEPS) {
          stepTime[plot_counter] = Time;
          WritePVD(outputFileName.c_str(), plot_counter);
        }

        FILE_LOGMatrixRM(DEBUGLOG, displacements, nNodes, ndim,
                         "Displacement Solution");
        writeFlag = false;
      }
      time_step_counter = time_step_counter + 1;
      dt = ExplicitTimeStepReduction * StableTimeStep();
      // Barrier not a must
      MPI_Barrier(MPI_COMM_WORLD);
    } // end explcit while loop
    FILE_LOG_MASTER(INFO, "End of Iterative Loop");

    // Write out the last time step
    // CustomPlot();
  } // end if ExplicitDynamic
  FILE_LOGMatrixRM(DEBUGLOG, displacements, nNodes, ndim,
                   "Final Displacement Solution");

  FinalizeFemTech();
  return 0;
}

void ApplyBoundaryConditions() {
  double tol = 1e-5;
  int count = 0;
  double AppliedDisp;

  AppliedDisp = Time * (dMax / tMax);

  // Apply Ramped Displacment
  // if (ExplicitDynamic || ImplicitDynamic) {
  //   if (Time < 0.1) {
  //     AppliedDisp = Time * (dMax / 0.1);
  //   } else {
  //     if (Time < 0.2) {
  //       AppliedDisp = dMax;
  //     } else {
  //       if (Time < 0.3) {
  //         AppliedDisp = dMax - (Time-0.2)*dMax/0.1;
  //       } else {
  //         if (Time < 0.4) {
  //           AppliedDisp = 0;
  //         } else {
  //           if (Time < 0.5) {
  //             AppliedDisp = (Time-0.4) * (dMax / 0.1);
  //           } else {
  //             if (Time < 0.6) {
  //               AppliedDisp = dMax;
  //             } else {
  //               if (Time < 0.7) {
  //                 AppliedDisp = dMax - (Time-0.6)*dMax/0.1;
  //               } else {
  //                 if (Time < 0.8) {
  //                   AppliedDisp = 0;
  //                 } else {
  //                   if (Time < 0.9) {
  //                     AppliedDisp = (Time-0.8) * (dMax / 0.1);
  //                   } else {
  //                     AppliedDisp = dMax;
  //                   }
  //                 }
  //               }
  //             }
  //           }
  //         }
  //       }
  //     }
  //   }
  // } else if (ImplicitStatic) {
  //   AppliedDisp = dMax;
  // }
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
      count = count + 1;
    }
    // if y coordinate = 0, constrain node to y plane (1-direction)
    index = ndim * i + 1;
    if (fabs(coordinates[index] - 0.0) < tol) {
      boundary[index] = 1;
      displacements[index] = 0.0;
      velocities[index] = 0.0;
      accelerations[index] = 0.0;
      count = count + 1;
    }
    // if z coordinate = 0, constrain node to z plane (2-direction)
    index = ndim * i + 2;
    if (fabs(coordinates[index] - 0.0) < tol) {
      boundary[index] = 1;
      displacements[index] = 0.0;
      velocities[index] = 0.0;
      accelerations[index] = 0.0;
      count = count + 1;
    }
    // if y coordinate = 1, apply disp. to node = 0.1 (1-direction)
    index = ndim * i + 1;
    if (fabs(coordinates[index] - 0.005) < tol) {
      boundary[index] = 1;
      count = count + 1;
      // note that this may have to be divided into
      // diplacement increments for both implicit and
      // explicit solver. In the future this would be
      // equal to some time dependent function i.e.,
      // CalculateDisplacement to get current increment out
      //  displacment to be applied.
      displacements[index] = AppliedDisp;
      velocities[index] = dMax / tMax;
      // For energy computations
      accelerations[index] = 0.0;
    }
  }
  FILE_LOG_MASTER(INFO, "Time = %10.5E, Applied Disp = %10.5E", Time,
                  AppliedDisp);
  return;
}

void ApplyVelocityBoundaryConditions(double) {
  double tol = 1e-5;
  int index;

  for (int i = 0; i < nNodes; i++) {
    // if x value = 0, constrain node to x plane (0-direction)
    index = ndim * i + 0;
    if (fabs(coordinates[index] - 0.0) < tol) {
      velocities_half[index] = 0.0;
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
    // if y coordinate = 1, apply disp. to node = 0.1 (1-direction)
    index = ndim * i + 1;
    if (fabs(coordinates[index] - 0.005) < tol) {
      velocities_half[index] = dMax / tMax;
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
    index = ndim * i + 0;
    if (fabs(coordinates[index] - 0.0) < tol) {
      boundary[index] = 1;
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
    // if y coordinate = 1, apply disp. to node = 0.1 (1-direction)
    index = ndim * i + 1;
    if (fabs(coordinates[index] - 0.005) < tol) {
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
      if (fabs(coordinates[ndim * i + x] - 0.005) < tol &&
          fabs(coordinates[ndim * i + y] - 0.005) < tol &&
          fabs(coordinates[ndim * i + z] - 0.005) < tol) {

        fprintf(datFile, "%13.5e  %13.5e  %13.5e  %13.5e  ", Time,
                displacements[ndim * i + x], displacements[ndim * i + y],
                displacements[ndim * i + z]);
      }
    }
    // Compute Cauchy stress for the element
    double stressElem[6];
    CalculateElementStress(0, stressElem);
    // Print all six stress components of 1st element stress
    fprintf(datFile, "%13.5e  %13.5e  %13.5e  %13.5e  %13.5e  %13.5e\n",
            stressElem[0], stressElem[1], stressElem[2], stressElem[5],
            stressElem[4], stressElem[3]);
  }
  fclose(datFile);
  return;
}
