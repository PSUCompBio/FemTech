#include "FemTech.h"
#include <math.h>
#include <string>

void ShellBMatrix(){
  for(int i=0; i<nshell; i++){
    Bshell[i*8+0] =  0.5*(hatcoordinates[connectivity[eptr[ShellID[i]]+1]*ndim+1] - hatcoordinates[connectivity[eptr[ShellID[i]]+3]*ndim+1])/areashell[i];
    Bshell[i*8+1] =  0.5*(hatcoordinates[connectivity[eptr[ShellID[i]]+2]*ndim+1] - hatcoordinates[connectivity[eptr[ShellID[i]]+0]*ndim+1])/areashell[i];
    Bshell[i*8+2] =  0.5*(-hatcoordinates[connectivity[eptr[ShellID[i]]+1]*ndim+1] + hatcoordinates[connectivity[eptr[ShellID[i]]+3]*ndim+1])/areashell[i];
    Bshell[i*8+3] =  0.5*(-hatcoordinates[connectivity[eptr[ShellID[i]]+2]*ndim+1] + hatcoordinates[connectivity[eptr[ShellID[i]]+0]*ndim+1])/areashell[i];
    Bshell[i*8+4] =  0.5*(-hatcoordinates[connectivity[eptr[ShellID[i]]+1]*ndim+0] + hatcoordinates[connectivity[eptr[ShellID[i]]+3]*ndim+0])/areashell[i];
    Bshell[i*8+5] =  0.5*(-hatcoordinates[connectivity[eptr[ShellID[i]]+2]*ndim+0] + hatcoordinates[connectivity[eptr[ShellID[i]]+0]*ndim+0])/areashell[i];
    Bshell[i*8+6] =  0.5*(hatcoordinates[connectivity[eptr[ShellID[i]]+1]*ndim+0] - hatcoordinates[connectivity[eptr[ShellID[i]]+3]*ndim+0])/areashell[i];
    Bshell[i*8+7] =  0.5*(hatcoordinates[connectivity[eptr[ShellID[i]]+2]*ndim+0] - hatcoordinates[connectivity[eptr[ShellID[i]]+0]*ndim+0])/areashell[i];
    hat_velocitystrain_half[i*3+0] = Bshell[i*8+0] * hat_velocities_half[connectivity[eptr[ShellID[i]]+0]*ndim] + Bshell[i*8+1] * hat_velocities_half[connectivity[eptr[ShellID[i]]+1]*ndim] + Bshell[i*8+2] * hat_velocities_half[connectivity[eptr[ShellID[i]]+2]*ndim] + Bshell[i*8+3] * hat_velocities_half[connectivity[eptr[ShellID[i]]+3]*ndim];
    hat_velocitystrain_half[i*3+1] = Bshell[i*8+4] * hat_velocities_half[connectivity[eptr[ShellID[i]]+0]*ndim+1] + Bshell[i*8+5] * hat_velocities_half[connectivity[eptr[ShellID[i]]+1]*ndim+1] + Bshell[i*8+6] * hat_velocities_half[connectivity[eptr[ShellID[i]]+2]*ndim+1] + Bshell[i*8+7] * hat_velocities_half[connectivity[eptr[ShellID[i]]+3]*ndim+1];
    hat_velocitystrain_half[i*3+2] = 0.5*(Bshell[i*8+4] * hat_velocities_half[connectivity[eptr[ShellID[i]]+0]*ndim] + Bshell[i*8+5] * hat_velocities_half[connectivity[eptr[ShellID[i]]+1]*ndim] + Bshell[i*8+6] * hat_velocities_half[connectivity[eptr[ShellID[i]]+2]*ndim] + Bshell[i*8+7] * hat_velocities_half[connectivity[eptr[ShellID[i]]+3]*ndim] + Bshell[i*8+0] * hat_velocities_half[connectivity[eptr[ShellID[i]]+0]*ndim+1] + Bshell[i*8+1] * hat_velocities_half[connectivity[eptr[ShellID[i]]+1]*ndim + 1] + Bshell[i*8+2] *  hat_velocities_half[connectivity[eptr[ShellID[i]]+2]*ndim + 1] + Bshell[i*8+3] * hat_velocities_half[connectivity[eptr[ShellID[i]]+3]*ndim + 1]);
  }
  return;
}
