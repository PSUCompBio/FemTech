#include "FemTech.h"

void CheckEnergy(void) {
  static double Wint_n = 0.0;
  static double Wext_n = 0.0;
  double sum = 0.0;
  for(int i = 0; i < nnodes*ndim; ++i) {
    // if(!boundary[i]) {
      sum += (displacements[i]-displacements_prev[i])*(fi[i]+fi_prev[i]);
    // }
  }
  sum *= 0.5;
  Wint_n += sum;
  printf("Internal Work : %.24f\n", Wint_n);
  double WKE = 0.0;
  for(int i = 0; i < nnodes*ndim; ++i) {
    // if(!boundary[i]) {
      WKE += mass[i]*velocities[i]*velocities[i];
    // }
  }
  WKE *= 0.5;
  printf("Kinetic Energy : %.24f\n", WKE);
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
    printf("\nERROR : Energy Violation\n\n");
  }
  // Wint_n = -WKE;
}
