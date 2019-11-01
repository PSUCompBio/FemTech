#include "FemTech.h"

/* Calculate center of mass of the whole body */
void GetBodyCenterofMass(double *cm) {
  for(int k = 0; k < ndim; ++k) {
    cm[k] = 0.0;
  }
  double totalMass = 0.0;
  double cmLocal[3];
  double volume, rho, massElem;
  int pide;
	for(int i = 0; i < nelements; ++i) {
    volume = CalculateCentroidAndVolume(i, cmLocal);
    pide = pid[i];
    rho = properties[MAXMATPARAMS * pide + 0];
    massElem = volume*rho;
    for (int j = 0; j < ndim; ++j) {
      cm[j] += massElem*cmLocal[j];
    }
    totalMass += massElem;
  }
  MPI_Allreduce(MPI_IN_PLACE, &totalMass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for (int j = 0; j < ndim; ++j) {
    MPI_Allreduce(MPI_IN_PLACE, &(cm[j]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    cm[j] = cm[j]/totalMass;
  }
}
