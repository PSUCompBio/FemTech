#include "FemTech.h"
#include "blas.h"

void CalculateEmbedDisp(int nodedof){
  double N1, N2, N3, N4, N5, N6, N7, N8, chi, eta, iota;
  chi = embedNC[nodedof*ndim+0];
  eta = embedNC[nodedof*ndim+1];
  iota = embedNC[nodedof*ndim+2];
  int embedid = int(nodedof/3);
  int host = embedinfo[embedid];
  N1 = ((1 - chi)*(1 - eta)*(1 - iota)) / 8;
  N2 = ((1 + chi)*(1 - eta)*(1 - iota)) / 8;
  N3 = ((1 + chi)*(1 + eta)*(1 - iota)) / 8;
  N4 = ((1 - chi)*(1 + eta)*(1 - iota)) / 8;
  N5 = ((1 - chi)*(1 - eta)*(1 + iota)) / 8;
  N6 = ((1 + chi)*(1 - eta)*(1 + iota)) / 8;
  N7 = ((1 + chi)*(1 + eta)*(1 + iota)) / 8;
  N8 = ((1 - chi)*(1 + eta)*(1 + iota)) / 8;
  int j = nodedof%3;
  displacements[nodedof] = N1*displacements[ndim*connectivity[eptr[host]]+j] + N2*displacements[ndim*connectivity[eptr[host]+1]+j] + N3*displacements[ndim*connectivity[eptr[host]+2]+j]
			 + N4*displacements[ndim*connectivity[eptr[host]+3]+j] + N5*displacements[ndim*connectivity[eptr[host]+4]+j] + N6*displacements[ndim*connectivity[eptr[host]+5]+j]
			 + N7*displacements[ndim*connectivity[eptr[host]+6]+j] + N8*displacements[ndim*connectivity[eptr[host]+7]+j];

if(nodedof==25){
//printf("%d %d %d %d\n", connectivity[eptr[host]], connectivity[eptr[host]+1] , connectivity[eptr[host]+2] , connectivity[eptr[host]+3]);
printf("%.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n",displacements[ndim*connectivity[host]+j], displacements[ndim*connectivity[host+1]+j], displacements[ndim*connectivity[host+2]+j], displacements[ndim*connectivity[host+3]+j], displacements[ndim*connectivity[host+4]+j], displacements[ndim*connectivity[host+5]+j], displacements[ndim*connectivity[host+6]+j], displacements[ndim*connectivity[host+7]+j]); 
}
  return;
}
