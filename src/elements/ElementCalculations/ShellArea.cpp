#include "FemTech.h"
#include <math.h>
#include <string>

void ShellArea(){
  for(int i = 0; i<nshell; i++){
    double elementCoordinates[12];
    for (int l = eptr[ShellID[i]], j = 0; l < eptr[ShellID[i]+1]; ++l, ++j) {
      for (int k = 0; k < 3; ++k) {
        int index = ndim*connectivity[l]+k;
        elementCoordinates[j*ndim+k] = coordinates[index]+displacements[index];
      }
    }
    const int index[4] = {0, 1, 2, 3};
    areashell[i] = areaHexahedronFace(elementCoordinates, &(index[0*4]));
//    areashell[i] = 0.5*((hatcoordinates[ndim*connectivity[eptr[ShellID[i]]+2]]-hatcoordinates[ndim*connectivity[eptr[ShellID[i]]]])*(hatcoordinates[ndim*connectivity[eptr[ShellID[i]]+3]+1]-hatcoordinates[ndim*connectivity[eptr[ShellID[i]]+1]+1])+(hatcoordinates[ndim*connectivity[eptr[ShellID[i]]+1]]-hatcoordinates[ndim*connectivity[eptr[ShellID[i]]+3]])*(hatcoordinates[ndim*connectivity[eptr[ShellID[i]]+2]+1]-hatcoordinates[ndim*connectivity[eptr[ShellID[i]]]+1]));
  }
  return;
}
