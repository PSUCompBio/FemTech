#include "FemTech.h"
#include <math.h>
#include <string>

void ShellHourglassing(){
  int h[4] = {1, -1, 1, -1};
  for(int i=0; i<nshell; i++){
    for(int j=0; j<4; j++){
    gamma_half[i*nshell+j] = h[j] - (h[0]*hatcoordinates[connectivity[eptr[ShellID[i]]+0]*ndim] + h[1]*hatcoordinates[connectivity[eptr[ShellID[i]]+1]*ndim] + h[2]*hatcoordinates[connectivity[eptr[ShellID[i]]+2]*ndim] + h[3]*hatcoordinates[connectivity[eptr[ShellID[i]]+3]*ndim+1])*Bshell[i*8+0*4+j] - (h[0]*hatcoordinates[connectivity[eptr[ShellID[i]]+0]*ndim+1] + h[1]*hatcoordinates[connectivity[eptr[ShellID[i]]+1]*ndim+1] + h[2]*hatcoordinates[connectivity[eptr[ShellID[i]]+2]*ndim+1] + h[3]*hatcoordinates[connectivity[eptr[ShellID[i]]+3]*ndim+1])*Bshell[i*8+1*4+j];
    }
    hg_strainrate[i*nshell + 0] = gamma_half[i*nshell+0]*hat_velocities_half[connectivity[eptr[ShellID[i]]+0]*ndim] + gamma_half[i*nshell+1]*hat_velocities_half[connectivity[eptr[ShellID[i]]+1]*ndim] + gamma_half[i*nshell+2]*hat_velocities_half[connectivity[eptr[ShellID[i]]+2]*ndim] + gamma_half[i*nshell+3]*hat_velocities_half[connectivity[eptr[ShellID[i]]+3]*ndim];
    hg_strainrate[i*nshell + 1] = gamma_half[i*nshell+0]*hat_velocities_half[connectivity[eptr[ShellID[i]]+0]*ndim+1] + gamma_half[i*nshell+1]*hat_velocities_half[connectivity[eptr[ShellID[i]]+1]*ndim+1] + gamma_half[i*nshell+2]*hat_velocities_half[connectivity[eptr[ShellID[i]]+2]*ndim+1] + gamma_half[i*nshell+3]*hat_velocities_half[connectivity[eptr[ShellID[i]]+3]*ndim+1];
  }
  return;
}
