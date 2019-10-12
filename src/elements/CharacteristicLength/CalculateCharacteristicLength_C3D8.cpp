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
#ifdef DEBUG
  if (debug && 1==0) {
    printf("Element %d volume : %12.6f\n", e, cl);
  }
#endif //DEBUG
  double areaMax = 0.0, faceArea;
  for (int i = 0; i < 6; ++i) {
    faceArea = areaHexahedronFace(elementCoordinates, &(index[i*4]));
#ifdef DEBUG
    if (debug && 1==0) {
      printf("Element %d, Side %d Area : %12.6f\n", e, i, faceArea);
    }
#endif //DEBUG
    if (faceArea > areaMax) {
      areaMax = faceArea;
    } 
  }
#ifdef DEBUG
  if (debug && 1==0) {
    printf("Element %d max Area : %12.6f\n", e, areaMax);
  }
#endif //DEBUG
  return cl/areaMax;
}
