#include "FemTech.h"
#include "blas.h"

/*Delare Functions*/
void ApplyBoundaryConditions();

double Time;
int nStep;
bool ImplicitStatic;
bool ImplicitDynamic;
bool ExplicitDynamic;

int main(int argc, char **argv){

	//true or false
  ImplicitStatic = false;
	ImplicitDynamic = false;
  ExplicitDynamic = true;

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

  if(debug) {
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
    ShapeFunctions();
    ReadMaterialProperties();
    ApplyBoundaryConditions();
    Assembly((char*)"stiffness");
    ApplySteadyBoundaryConditions();
    SolveSteadyImplicit();
    Time = 1.0;
    nStep = 1;
    /* Write final, deformed configuration*/
    WriteVTU(argv[1], nStep, Time);
  }
	else if (ImplicitDynamic) {
    // Dynamic Implicit solution using Newmark's scheme for time integration
    double dt = 0.1;
    double tMax = 10.0;
    ShapeFunctions();
    ReadMaterialProperties();
    ApplyBoundaryConditions();
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
    double tMax = 1e-3; //in seconds
		double t=0.0;
    ShapeFunctions();
    ReadMaterialProperties();
    ApplyBoundaryConditions();
    Assembly((char*)"stiffness");
    /*  Step-1: Calculate the mass matrix similar to that of belytschko. */
    Assembly((char*)"mass"); // Add Direct-lumped as an option
     /* Step-2: getforce step from Belytschko */
	//	getforce()
	//	dt = fe_getTimeStep();
		/* Step-3: Calculate accelerations */
  //  fe_calculateAccln(accelerations, m_system, F_net);
   // displacements_prev = displacements;
		/* Step-4: Time loop starts....*/
	//	time_step_counter = time_step_counter + 1;
	//	clock_t s, s_prev, ds;
	//	s = clock();
		while (t <= tMax) {
 			t=t+dt;
			printf("t = %3.3e\n",t);

		/** Steps - 4,5,6 and 7 from Belytschko Box 6.1 - Update time, velocity and displacements */
		//fe_timeUpdate(U, V, V_half, A, t, dT, "newmark-beta-central-difference");

	  /** Update Loading Conditions - time dependent loading conditions */
	  //fe_apply_bc_load(fe, t);

		/** Step - 8 from Belytschko Box 6.1 - Calculate net nodal force*/
    //fe_getforce(F_net, ndof, U, fe, time_step_counter, U_prev, dT, f_damp_curr, d_static, d_fatigue, d_tot, lambda_min, lambda_max, lambda_min_cycle, lambda_max_cycle, d_avg, n_load_cycle_full, n_load_cycle_partial, t, t_plot); // Calculating the force term.

    /** Step - 9 from Belytschko Box 6.1 - Calculate Accelerations */
    //fe_calculateAccln(A, m_system, F_net); // Calculating the new accelerations from total nodal forces.
    //fe_apply_bc_acceleration(A, t);

    /** Step- 10 from Belytschko Box 6.1 - Second Partial Update of Nodal Velocities */
    //fe_timeUpdate_velocity(V, V_half, A, t, dT, "newmark-beta-central-difference");

    //fi_curr = fe - F_net;
    //fe_calculateFR(fr_curr, fi_curr, m_system, A);

    /** Step - 11 from Belytschko Box 6.1 - Calculating energies and Checking Energy Balance */
    //fe_checkEnergies(U_prev, U, fi_prev, fi_curr, f_damp_prev, f_damp_curr, fe_prev, fe, fr_prev, fr_curr, m_system, V, energy_int_old, energy_int_new, energy_vd_old, energy_vd_new, energy_ext_old, energy_ext_new, energy_k

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

void ApplyBoundaryConditions(){
  double tol = 1e-5;
  int count = 0;
  printf("DEBUG : \n");

  for(int i=0;i<nnodes;i++){
    // if x value = 0, constrain node to x plane (0-direction)
    if(fabs(coordinates[ndim*i+0]-0.0) <tol){
      boundary[ndim*i+0]=1;
      printf("node : %d x : %d\n", i, count);
      count = count+1;
    }
    // if y coordinate = 0, constrain node to y plane (1-direction)
    if(fabs(coordinates[ndim*i+1]-0.0) <tol){
      boundary[ndim*i+1]=1;
      printf("node : %d y : %d\n", i, count);
      count = count+1;
    }
    // if z coordinate = 0, constrain node to z plane (2-direction)
    if(fabs(coordinates[ndim*i+2]-0.0) <tol){
      boundary[ndim*i+2]=1;
      printf("node : %d z : %d\n", i, count);
      count = count+1;
    }
    // if y coordinate = 1, apply disp. to node = 0.1 (1-direction)
    if(fabs(coordinates[ndim*i+1]-1.0) <tol){
      boundary[ndim*i+1]=1;
      printf("node : %d y2 : %d\n", i, count);
      count = count+1;
      // note that this may have to be divided into
      // diplacement increments for both implicit and
      // explicit solver. In the future this would be
      // equal to some time dependent function i.e.,
      // CalculateDisplacement to get current increment out
      //  displacment to be applied.
      displacements[ndim*i+1] = 0.1;
    }
  }
  printf("DEBUG : \n");
  return;
}
