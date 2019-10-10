#include "FemTech.h"

void CheckEnergy(double time) {
  static double Wint_n = 0.0;
  static double Wext_n = 0.0;

  /* Energy Calculations are based on Belytschko et al. book
		 Section 6.2.3 Energy Balance */

  /* Calculate Internal and External Energies */
  /* Calculate Kinetic Energy */
  double sum_Wint_n = 0.0;
	double sum_Wext_n = 0.0;
	double delta_d = 0.0;
  double WKE = 0.0;

  // Calculate in parallel
  // Loop over all the nodes
  for(int i = 0; i < nnodes; ++i) {
    bool includeSum = true;
    for (int j = 0; (j < sendProcessCount) && includeSum; ++j) {
      // If sending to a node with lower rank its already included
      // in the KE calculation
      if (sendProcessID[j] < world_rank) {
        for (int k = sendNeighbourCountCum[j]; k < sendNeighbourCountCum[j+1]; ++k) {
          if (sendNodeIndex[k] == i) {
            includeSum = false;
            break;
            // skip addition of KE from this node
          }
        }
      }
    }
    if (includeSum) {
      int index = i*ndim;
      for (int j = 0; j < ndim; ++j) {
        int indexJ = index + j;
			  delta_d = displacements[indexJ]-displacements_prev[indexJ];
        WKE += mass[indexJ]*velocities[indexJ]*velocities[indexJ];
        if (boundary[indexJ]) {
          // Fext = Fint + m*Acceleration
          double reaction = fi_prev[indexJ] + fi[indexJ] + mass[indexJ]*(accelerations[indexJ]+accelerations_prev[indexJ]);

			    sum_Wext_n += delta_d*reaction;
        }
        sum_Wint_n += delta_d*(fi_prev[indexJ] + fi[indexJ]); /* equ 6.2.14 */
			  sum_Wext_n += delta_d*(fe_prev[indexJ] + fe[indexJ]); /* equ 6.2.15 */
      }
    }
  }
  WKE *= 0.5;
  sum_Wint_n *= 0.5;
  sum_Wext_n *= 0.5;

  double WKE_Total = 0.0, Wint_n_total = 0.0, Wext_n_total = 0.0;
  MPI_Reduce(&WKE, &WKE_Total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sum_Wint_n, &Wint_n_total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sum_Wext_n, &Wext_n_total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (world_rank ==0) {
    Wint_n += Wint_n_total;
    Wext_n += Wext_n_total;
    double total = fabs(WKE_Total+Wint_n-Wext_n);
#ifdef DEBUG
	if(debug && 1 == 0){
  	printf("Internal Work : %15.9e\n", Wint_n);
  	printf("External Work : %15.9e\n", Wext_n);
    printf("Kinetic Energy : %15.9e\n", WKE_Total);
    printf("Total Energy : %15.9e\n", total);
}
#endif //DEBUG

    double max = fabs(Wint_n);
    if (max < fabs(Wext_n)) {
      max = fabs(Wext_n);
    }
    if (max < fabs(WKE_Total)) {
      max = fabs(WKE_Total);
    }
    const double epsilon = 0.01;
    if (total > epsilon*max) {
      printf("\nERROR - Energy Violation:  Total = %15.9e, Max = %15.9e, Error%% : %10.2f \n", total, max, total*100.0/max);
    }
    fprintf(energyFile, "%12.6e %12.6e  %12.6e  %12.6e %12.6e\n", time,
            Wint_n, Wext_n, WKE_Total, total);
  }
}
