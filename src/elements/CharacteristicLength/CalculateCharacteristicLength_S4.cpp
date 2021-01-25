#include "FemTech.h"

double CalculateCharacteristicLength_S4(int e) {
  double elementCoordinates[12];
  for (int i = eptr[e], j = 0; i < eptr[e+1]; ++i, ++j) {
    for (int k = 0; k < 3; ++k) {
      int index = ndim*connectivity[i]+k;
      elementCoordinates[j*ndim+k] = coordinates[index]+displacements[index];
    }
  }
  const int index[4] = {0, 1, 2, 3};

  double faceArea;
  faceArea = areaHexahedronFace(elementCoordinates, &(index[0*4]));
  FILE_LOG_SINGLE(DEBUGLOGIGNORE, "Element %d, Area : %12.6f", e, faceArea);

 const int indexlength[8] = {0, 1, 1, 2, 2, 3, 3, 0};

  double maxlength = 0.0, length;
  for(int i=0; i<4; i++){
    length = shellsidelength(elementCoordinates, &(index[i*2]));
    if (length > maxlength) {
    maxlength = length;
   }
 }

  FILE_LOG_SINGLE(DEBUGLOGIGNORE, "Element %d max length : %12.6f", e, maxlength);
  return faceArea/maxlength;
}
