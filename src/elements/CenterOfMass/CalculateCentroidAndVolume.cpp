#include "FemTech.h"

double CalculateCentroidAndVolume(int e, double *cm) {
  double volume;
  if (strcmp(ElementType[e], "C3D8") == 0) {
    volume = CalculateCentroidAndVolume_C3D8(e, cm);
  } else {
    if (strcmp(ElementType[e], "C3D4") == 0) {
      CalculateCentroid_C3D4(e, cm);
      double coord[12];
      for (int i = eptr[e], j = 0; i < eptr[e+1]; ++i, ++j) {
        for (int k = 0; k < 3; ++k) {
          int index = ndim*connectivity[i]+k;
          coord[j*ndim+k] = coordinates[index]+displacements[index];
        }
      }
      volume = volumeTetrahedron(coord);
    } else {
      FILE_LOG_SINGLE(ERROR, "Unknown Element Type Encountered");
      TerminateFemTech(3);
    }
  }
  return volume;
}
