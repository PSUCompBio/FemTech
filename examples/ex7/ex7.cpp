#include "FemTech.h"
#include "blas.h"

/*Delare Functions*/
void ApplyBoundaryConditions(double Time,double dMax, double tMax);

/* Global Variables/Parameters  - could be moved to parameters.h file?  */
double Time;
int nStep;
int nSteps;
int nPlotSteps = 10;
bool ImplicitStatic = false;
bool ImplicitDynamic = true;
bool ExplicitDynamic = false;
double ExplicitTimeStepReduction = 0.8;
double FailureTimeStep = 1e-11;

int main(int argc, char **argv){

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

  if(debug && 1==0) {
    // Printing local arrays of processor (this section can be removed)
    printf("\neptr array in processor %d after partitioning = ", world_rank);
    for (int i = 0; i <= nelements; i++) {
      printf("%d ", eptr[i]);
    }
    printf("\n");
    printf("\neind array in processor %d after partitioning =", world_rank);
    for (int i = 0; i < nelements; i++) {
      printf(" (%d)  ", i);
      for (int j = eptr[i]; j < eptr[i + 1]; j++) {
        printf("%d ", connectivity[j]);
      }
    }
    printf("\n");
    printf("\nType/PartID of element in processor %d after partitioning = ", world_rank);
    for (int i = 0; i < nelements; i++) {
      printf("%s/%d  ", ElementType[i], pid[i]);
    }
    printf("\n");
    printf("\nSize of coordinates array in processor %d after partitioning = %d\n", world_rank, nnodes * ndim);
    printf("\nCoordinates array in processor %d after partitioning =", world_rank);
    for (int i = 0; i < nnodes; i++) {
      printf(" (%d)  ", i);
      for (int j = 0; j < ndim; j++) {
        printf("%.*f ", 1, coordinates[ndim * i + j]);
      }
    }
    printf("\n");
  }

  /* Write inital, undeformed configuration*/
	Time = 0.0;
	nStep = 0;
  WriteVTU(argv[1], nStep, Time);

  if (ImplicitStatic) {
    // Static solution
    double dMax = 0.05; // max displacment in meters
		double tMax = 1.0;
    ShapeFunctions();
    ReadMaterialProperties();
    Time = 1.0;
    ApplyBoundaryConditions(Time,dMax,tMax);
    Assembly((char*)"stiffness");
    ApplySteadyBoundaryConditions();
    SolveSteadyImplicit();
    Time = tMax;
    nStep = 1;
    /* Write final, deformed configuration*/
    WriteVTU(argv[1], nStep, Time);
    if (debug) {
      printf("DEBUG : Printing Displacement Solution\n");
      for (int i = 0; i < nnodes; ++i) {
        for (int j = 0; j < ndim; ++j) {
          printf("%12.4f", displacements[i*ndim+j]);
        }
        printf("\n");
      }
  }

  }
	else if (ImplicitDynamic) {
    // Dynamic Implicit solution using Newmark's scheme for time integration
    double dt = 0.05;
    double tMax = 10.0;
    double dMax = 0.05;// max displacment in meters
    ShapeFunctions();
    ReadMaterialProperties();
    Time = tMax;
    ApplyBoundaryConditions(Time,dMax,tMax);
    Time = 0.0;
    Assembly((char*)"stiffness");
    Assembly((char*)"mass");
    /* beta and gamma of Newmark's scheme */
    double beta = 0.25;
    double gamma = 0.5;
    SolveUnsteadyNewmarkImplicit(beta, gamma, dt, tMax, argv[1]);
  }
	else if (ExplicitDynamic){
		// Dynamic Explcit solution using....
		double dt = 1e-4;
    double tMax = 1e-3; //max simulation time in seconds
    double dMax = 0.1; // max displacment in meters
		double Time=0.0;
		int time_step_counter = 0;
		int plot_counter = 0;
		/** Central Difference Method - Beta and Gamma */
		double beta = 0;
		double gamma = 0.5;

    ShapeFunctions();
    ReadMaterialProperties();
    //ApplyBoundaryConditions(Time,dMax,tMax);
    Assembly((char*)"stiffness");
    /*  Step-1: Calculate the mass matrix similar to that of belytschko. */

    Assembly((char*)"mass"); // Add Direct-lumped as an option
		/* Step-2: getforce step from Belytschko */

    GetForce(); // Calculating the force term.

    /* obtain stable time step */
   	dt = ExplicitTimeStepReduction*StableTimeStep();
		nSteps = (int)(tMax/dt);
		int nsteps_plot = (int)(nSteps/nPlotSteps);
		printf("inital dt = %3.3e, nSteps = %d, nsteps_plot = %d\n",dt,nSteps,nsteps_plot);

		/* Step-3: Calculate accelerations */
   	CalculateAccelerations();

 		// Save old displacments
		for(int i=0;i<ndim*nnodes;i++){
   		displacements_prev[i] = displacements[i];
		}

		/* Step-4: Time loop starts....*/
		time_step_counter = time_step_counter + 1;
	//	clock_t s, s_prev, ds;
	//	s = clock();
		//while (Time <= tMax) {
		for(int k=0;k<1;k++){
 			Time=Time+dt; /*Update the time by adding full time step */

			/* Steps - 4,5,6 and 7 from Belytschko Box 6.1 - Update time, velocity and displacements */
			/* Partially Update Nodal Velocities */
			for(int i=0;i<ndim*nnodes;i++){
				velocities_half[i] = velocities[i] + ((1.0 - gamma)*dt*accelerations[i]);
			}

			/* Update Nodal Displacements */
			for(int i=0;i<ndim*nnodes;i++){
				displacements[i] = displacements[i]+
													(dt*velocities[i]) +
													((pow(dt,2)/2.0)*(1.0-(2.0*beta))*accelerations[i]) +
													(beta*dt*dt*accelerations[i]);
			}
	 		/** Update Loading Conditions - time dependent loading conditions */
			ApplyBoundaryConditions(Time,dMax,tMax);

			/* Step - 8 from Belytschko Box 6.1 - Calculate net nodal force*/
    	GetForce(); // Calculating the force term.

    	/* Step - 9 from Belytschko Box 6.1 - Calculate Accelerations */
    	CalculateAccelerations(); // Calculating the new accelerations from total nodal forces.
    	//fe_apply_bc_acceleration(A, t);

    	/** Step- 10 from Belytschko Box 6.1 - Second Partial Update of Nodal Velocities */
    	// using newmark-beta-central-difference
			for(int i=0;i<ndim*nnodes;i++){
				velocities[i] = velocities_half[i] + (gamma*dt*accelerations[i]);
			}
			//fe_apply_bc_velocity(V, t);

    //fi_curr = fe - F_net;
			for(int i=0;i<ndim*nnodes;i++){
				fi_curr[i] = fe[i] - f_net[i];
			}

      CalculateFR();

    /** Step - 11 from Belytschko Box 6.1 - Calculating energies and Checking Energy Balance */
    //fe_checkEnergies(U_prev, U, fi_prev, fi_curr, f_damp_prev, f_damp_curr, fe_prev, fe, fr_prev, fr_curr, m_system, V, energy_int_old, energy_int_new, energy_vd_old, energy_vd_new, energy_ext_old, energy_ext_new, energy_k
    //printf("mod: %d\n", time_step_counter % nsteps_plot);
		if(time_step_counter % nsteps_plot == 0 ){
			plot_counter = plot_counter + 1;

			printf("------Plot %d: WriteVTU\n",plot_counter);
			WriteVTU(argv[1], plot_counter, Time);
		}
		time_step_counter = time_step_counter + 1;
		dt = StableTimeStep();

		} // end explcit while loop
	}// end if ExplicitDynamic

	/* Below are things to do at end of program */
	if(world_rank == 0){
		WritePVD(argv[1], nStep, Time);
	}
  FreeArrays();
  MPI_Finalize();
  return 0;
}

void ApplyBoundaryConditions(double Time, double dMax, double tMax){
  double tol = 1e-5;
  int count = 0;
	double AppliedDisp;

	// Apply Ramped Displacment
	if (ExplicitDynamic || ImplicitDynamic){
		AppliedDisp = Time*(dMax/tMax);
  }
  else if (ImplicitStatic){
    AppliedDisp = dMax;
  }


  for(int i=0;i<nnodes;i++){
    // if x value = 0, constrain node to x plane (0-direction)
    if(fabs(coordinates[ndim*i+0]-0.0) <tol){
      boundary[ndim*i+0]=1;
      //printf("node : %d x : %d\n", i, count);
      count = count+1;
    }
    // if y coordinate = 0, constrain node to y plane (1-direction)
    if(fabs(coordinates[ndim*i+1]-0.0) <tol){
      boundary[ndim*i+1]=1;
      //printf("node : %d y : %d\n", i, count);
      count = count+1;
    }
    // if z coordinate = 0, constrain node to z plane (2-direction)
    if(fabs(coordinates[ndim*i+2]-0.0) <tol){
      boundary[ndim*i+2]=1;
      //printf("node : %d z : %d\n", i, count);
      count = count+1;
    }
    // if y coordinate = 1, apply disp. to node = 0.1 (1-direction)
    if(fabs(coordinates[ndim*i+1]-1.0) <tol){
      boundary[ndim*i+1]=1;
      //printf("node : %d y2 : %d\n", i, count);
      count = count+1;
      // note that this may have to be divided into
      // diplacement increments for both implicit and
      // explicit solver. In the future this would be
      // equal to some time dependent function i.e.,
      // CalculateDisplacement to get current increment out
      //  displacment to be applied.
      displacements[ndim*i+1] = AppliedDisp;
    }
  }
  //printf("Time = %3.3e, Applied Disp = %3.3e\n",Time,AppliedDisp);
  return;
}
