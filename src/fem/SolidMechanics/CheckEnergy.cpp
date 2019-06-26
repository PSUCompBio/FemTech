#include "FemTech.h"

void CheckEnergy(void) {
  static double Wint_n = 0.0;
  static double Wext_n = 0.0;

  /* Energy Calculations are based on Belytschko et al. book
		 Section 6.2.3 Energy Balance */

  /* Calculate Internal and External Energies */
  double sum_Wint_n = 0.0;
	double sum_Wext_n = 0.0;
	double delta_d=0.0;

  for(int i = 0; i < nnodes*ndim; ++i) {
			delta_d=displacements[i]-displacements_prev[i];
      sum_Wint_n += delta_d*(fi_prev[i] + fi[i]); /* equ 6.2.14 */
			sum_Wext_n += delta_d*(fe_prev[i] + fe[i]); /* equ 6.2.15 */
  }

  sum_Wint_n *= 0.5;
  sum_Wext_n *= 0.5;
  Wint_n += sum_Wint_n;
  Wext_n += sum_Wext_n;

#ifdef DEBUG
	if(debug){
  	printf("Internal Work : %.24f\n", Wint_n);
	}
#endif //DEBUG
  
/* Calculate Kinetic Energy */
  double WKE = 0.0;
  for(int i = 0; i < nnodes*ndim; ++i) {
      WKE += mass[i]*velocities[i]*velocities[i];
  }
  WKE *= 0.5;
  
#ifdef DEBUG
	if(debug){
  	printf("Kinetic Energy : %.24f\n", WKE);
	}
#endif //DEBUG
  
  double total = fabs(WKE+Wint_n-Wext_n);

  double max = fabs(Wint_n);
  if (max < fabs(Wext_n)) {
    max = fabs(Wext_n);
  }
  if (max < fabs(WKE)) {
    max = fabs(WKE);
  }
  const double epsilon = 0.01;
  if (total > epsilon*max) {
    printf("\nERROR - Energy Violation:  IW = %3.3e, KE=%3.3e \n\n",Wint_n,WKE);
  }
  // Wint_n = -WKE;
}
