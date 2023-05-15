#include "FemTech.h"

double CalculateVolume_T3D2(int e) {
  double volume = 0.0;
  double length;
  double area = 1.0; //temporary
  double elementCoordinates[6];
  for (int i = eptr[e], j = 0; i < eptr[e+1]; ++i, ++j) {
    for (int k = 0; k < 3; ++k) {
      int index = ndim*connectivity[i]+k;
      elementCoordinates[j*ndim+k] = coordinates[index]+displacements[index];
    }
  }
  length = std::abs(elementCoordinates[1]-elementCoordinates[0]);
  volume = area*length;
  return volume;
}
