#include "FemTech.h"

double distancePointPlane(double *x0, double* x1, double* x2, double* x3);

double CalculateCharacteristicLength_C3D4(int e) {
  int index[16] = {0, 1, 2, 3, 1, 2, 3, 0, 2, 3, 0, 1, 3, 0, 1, 2};
  double elementCoordinates[12];
  for (int i = eptr[e], j = 0; i < eptr[e+1]; ++i, ++j) {
    for (int k = 0; k < 3; ++k) {
      int index = ndim*connectivity[i]+k;
      elementCoordinates[j*ndim+k] = coordinates[index]+displacements[index];
    }
  }
  double altitudeMin = 1e6, minAlt;
  for (int i = 0; i < 4; ++i) {
    minAlt = distancePointPlane(&(elementCoordinates[index[4*i]*ndim]), \
        &(elementCoordinates[index[4*i+1]*ndim]), \
        &(elementCoordinates[index[4*i+2]*ndim]), \
        &(elementCoordinates[index[4*i+3]*ndim])); 
    if (debug) {
      printf("Element %d, altitude %d : %12.6f\n", e, i, minAlt);
    }
    if (minAlt < altitudeMin) {
      altitudeMin = minAlt;
    }
  }
  if (debug) {
    printf("Element %d, min altitude : %12.6f\n", e, altitudeMin);
  }
  return altitudeMin;
}

double distancePointPlane(double *x0, double* x1, double* x2, double* x3) {
  // Computed distance of point x0 from plane containing x1, x2, x3
  double v1[3], v2[3], normal[3];
  for (int i = 0; i < 3; ++i) {
    v1[i] = x2[i]-x1[i];
    v2[i] = x3[i]-x1[i];
  }
  crossProduct(v1, v2, normal);
  double crossNorm = norm3D(normal);
  for (int i = 0; i < 3; ++i) {
    normal[i] /= crossNorm;
  }
  double v3[3];
  for (int i = 0; i < 3; ++i) {
    v3[i] = x0[i] - x1[i];
  }
  return fabs(dotProduct3D(normal, v3));
}
