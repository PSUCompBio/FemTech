#include "FemTech.h"
#include "blas.h"

double CalculateStrain(int e){
  // computing Cavg
  double *Cavg = (double *)malloc(ndim * ndim * sizeof(double)); //right cauchy green deformation tensor avg
  dgemm_(chy, chn, &ndim, &ndim, &ndim, &one, Favg, &ndim, Favg, &ndim, &zero, Cavg, &ndim); //Cavg = F^T * F
  //making the identity matrix
  double *I = (double *)malloc(ndim * ndim * sizeof(double));
  for (int i = 0; i < ndim*ndim; i++){ //maybe there is a better way to make this matrix
      I[i] = 0.0;
      }
  for (int i = 0; i < ndim; i++){
      I[i*ndim + i] = 1.0;
      }
  //Computing Strain
  double alpha = -1.0;
  int New = ndim*ndim;
  daxpy_(&New, &alpha, I, &oneI, Cavg, &oneI); // Cavg = Cavg - I (this is an intermediate step to calculate strains, has no physical meaning)
  double max = -1;
  for(int i = 0; i < ndim*ndim; i++){
    Eavg[e*ndim*ndim + i] = 0.5*(Cavg[i]); // Eavg = 0.5 * (Cavg - I)
    if (max < fabs(0.5*(Cavg[i]))) max = fabs(0.5*(Cavg[i]));
  }

  free(Cavg);
  free(I);
  return max;
}
