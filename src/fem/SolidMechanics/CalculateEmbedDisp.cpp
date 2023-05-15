#include "FemTech.h"
#include "blas.h"

void CalculateEmbedDisp(int embedid, int dirn){
  double N1, N2, N3, N4, N5, N6, N7, N8, chi, eta, iota;
  chi = embedNC[embedid*ndim+0];
  eta = embedNC[embedid*ndim+1];
  iota = embedNC[embedid*ndim+2];
  int host = nodeconstrain[embedid];
  N1 = ((1 - chi)*(1 - eta)*(1 - iota)) / 8;
  N2 = ((1 + chi)*(1 - eta)*(1 - iota)) / 8;
  N3 = ((1 + chi)*(1 + eta)*(1 - iota)) / 8;
  N4 = ((1 - chi)*(1 + eta)*(1 - iota)) / 8;
  N5 = ((1 - chi)*(1 - eta)*(1 + iota)) / 8;
  N6 = ((1 + chi)*(1 - eta)*(1 + iota)) / 8;
  N7 = ((1 + chi)*(1 + eta)*(1 + iota)) / 8;
  N8 = ((1 - chi)*(1 + eta)*(1 + iota)) / 8;
  displacements[ndim*embedid+dirn] = N1*displacements[ndim*connectivity[eptr[host]]+dirn] + N2*displacements[ndim*connectivity[eptr[host]+1]+dirn] + N3*displacements[ndim*connectivity[eptr[host]+2]+dirn]
			 + N4*displacements[ndim*connectivity[eptr[host]+3]+dirn] + N5*displacements[ndim*connectivity[eptr[host]+4]+dirn] + N6*displacements[ndim*connectivity[eptr[host]+5]+dirn]
			 + N7*displacements[ndim*connectivity[eptr[host]+6]+dirn] + N8*displacements[ndim*connectivity[eptr[host]+7]+dirn];
  accelerations[ndim*embedid+dirn] = N1*accelerations[ndim*connectivity[eptr[host]]+dirn] + N2*accelerations[ndim*connectivity[eptr[host]+1]+dirn] + N3*accelerations[ndim*connectivity[eptr[host]+2]+dirn]
			 + N4*accelerations[ndim*connectivity[eptr[host]+3]+dirn] + N5*accelerations[ndim*connectivity[eptr[host]+4]+dirn] + N6*accelerations[ndim*connectivity[eptr[host]+5]+dirn]
			 + N7*accelerations[ndim*connectivity[eptr[host]+6]+dirn] + N8*accelerations[ndim*connectivity[eptr[host]+7]+dirn];
  return;
}
