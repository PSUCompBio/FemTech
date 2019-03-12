#include "FemTech.h"

void GetForce(double timeFinal) {
  double tCritical = timeFinal;
}

void ExplicitDynamics(double timeFinal, char* name) {
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
		while (t <= timeFinal) {
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
}
