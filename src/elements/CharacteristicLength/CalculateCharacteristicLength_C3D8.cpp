#include "FemTech.h"

double CalculateCharacteristicLength_C3D8(int e) {
  double cl;
  double elementCoordinates[24];
  for (int i = eptr[e], j = 0; i < eptr[e+1]; ++i, ++j) {
    for (int k = 0; k < 3; ++k) {
      int index = ndim*connectivity[i]+k;
      elementCoordinates[j*ndim+k] = coordinates[index]+displacements[index];
    }
  }
  const int index[24] = {0, 1, 2, 3, 4, 5, 6, 7, 0, 3, 7, 4, 1, 2, 6, 5, \
                   0, 1, 5, 4, 3, 2, 6, 7};
  cl = volumeHexahedron(elementCoordinates);

  FILE_LOG_SINGLE(DEBUGLOGIGNORE, "Element %d volume : %12.6f", e, cl);

  double areaMax = 0.0, faceArea;
  for (int i = 0; i < 6; ++i) {
    faceArea = areaHexahedronFace(elementCoordinates, &(index[i*4]));

    FILE_LOG_SINGLE(DEBUGLOGIGNORE, "Element %d, Side %d Area : %12.6f", e, i, faceArea);

    if (faceArea > areaMax) {
      areaMax = faceArea;
    } 
  }

  FILE_LOG_SINGLE(DEBUGLOGIGNORE, "Element %d max Area : %12.6f", e, areaMax);

  return cl/areaMax;
}
