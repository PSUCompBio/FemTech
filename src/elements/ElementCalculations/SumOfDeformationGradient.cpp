#include "FemTech.h"
#include "blas.h"

void SumOfDeformationGradient(int e, int gp){
  int index = fptr[e] + ndim * ndim * gp;
  double *F_gp = &(F[index]);
  for(int i = 0; i < ndim*ndim; i++)
     {Favg[e*ndim*ndim + i] += F_gp[i];} //adding F components of all gauss points
}
