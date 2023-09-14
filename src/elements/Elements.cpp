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

/* Calculate volume of each partID */
void computePartVolume(double *volume, double *elemVolume) {
  int pide;
  double volumeE;
  for (int i = 0; i < nelements; ++i) {
    if(pid[i]==10){
	elemVolume[i] = 0.0;
	volume[10] = 0.0;}
    else{
    volumeE = calculateVolume(i);
    pide = pid[i];
    volume[pide] = volume[pide] + volumeE;
    elemVolume[i] = volumeE;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, volume, nPIDglobal, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}
