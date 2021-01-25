#include "FemTech.h"
#include <math.h>
#include <string>

void CorotationalSystemTransform(){

  for(int i=0; i<nshell; i++){
    for(int j=0; j<4; j++){
      hatcoordinates[connectivity[eptr[ShellID[i]]+j]*ndim] = corotationalx[i*ndim] * (coordinates[connectivity[eptr[ShellID[i]]+j]*ndim]+displacements[connectivity[eptr[ShellID[i]]+j]*ndim]) + corotationalx[i*ndim+1] * (coordinates[connectivity[eptr[ShellID[i]]+j]*ndim+1]+displacements[connectivity[eptr[ShellID[i]]+j]*ndim+1]) + corotationalx[i*ndim+2] * (coordinates[connectivity[eptr[ShellID[i]]+j]*ndim+2]+displacements[connectivity[eptr[ShellID[i]]+j]*ndim+2]);
      hatcoordinates[connectivity[eptr[ShellID[i]]+j]*ndim+1] = corotationaly[i*ndim] * (coordinates[connectivity[eptr[ShellID[i]]+j]*ndim]+displacements[connectivity[eptr[ShellID[i]]+j]*ndim]) + corotationaly[i*ndim+1] * (coordinates[connectivity[eptr[ShellID[i]]+j]*ndim+1]+displacements[connectivity[eptr[ShellID[i]]+j]*ndim+1]) + corotationaly[i*ndim+2] * (coordinates[connectivity[eptr[ShellID[i]]+j]*ndim+2]+displacements[connectivity[eptr[ShellID[i]]+j]*ndim+2]);
      hatcoordinates[connectivity[eptr[ShellID[i]]+j]*ndim+2] = corotationalz[i*ndim] * (coordinates[connectivity[eptr[ShellID[i]]+j]*ndim]+displacements[connectivity[eptr[ShellID[i]]+j]*ndim]) + corotationalz[i*ndim+1] * (coordinates[connectivity[eptr[ShellID[i]]+j]*ndim+1]+displacements[connectivity[eptr[ShellID[i]]+j]*ndim+1]) + corotationalz[i*ndim+2] * (coordinates[connectivity[eptr[ShellID[i]]+j]*ndim+2]+displacements[connectivity[eptr[ShellID[i]]+j]*ndim+2]);
      hat_velocities_half[connectivity[eptr[ShellID[i]]+j]*ndim] = corotationalx[i*ndim] * velocities_half[connectivity[eptr[ShellID[i]]+j]*ndim] + corotationalx[i*ndim+1] * velocities_half[connectivity[eptr[ShellID[i]]+j]*ndim+1] + corotationalx[i*ndim+2] * velocities_half[connectivity[eptr[ShellID[i]]+j]*ndim+2];
      hat_velocities_half[connectivity[eptr[ShellID[i]]+j]*ndim+1] = corotationaly[i*ndim] * velocities_half[connectivity[eptr[ShellID[i]]+j]*ndim] + corotationaly[i*ndim+1] * velocities_half[connectivity[eptr[ShellID[i]]+j]*ndim+1] + corotationaly[i*ndim+2] * velocities_half[connectivity[eptr[ShellID[i]]+j]*ndim+2];;
      hat_velocities_half[connectivity[eptr[ShellID[i]]+j]*ndim+2] = corotationalz[i*ndim] * velocities_half[connectivity[eptr[ShellID[i]]+j]*ndim] + corotationalz[i*ndim+1] * velocities_half[connectivity[eptr[ShellID[i]]+j]*ndim+1] + corotationalz[i*ndim+2] * velocities_half[connectivity[eptr[ShellID[i]]+j]*ndim+2];
    }
  }

  return;
}
