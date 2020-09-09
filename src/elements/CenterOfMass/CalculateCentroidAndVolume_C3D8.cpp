#include "FemTech.h"

double CalculateCentroidAndVolume_C3D8(int e, double *cm) {
  double volume = 0.0, volumeTet;
  for (int k = 0; k < ndim; ++k) {
    cm[k] = 0.0;
  }
  // Split hexahedron to 6 tetrahedrons and compute centroid and volume
  // Each set of 4 indices forms a tetrahedron from hexahedron coords
  const int index[24] = {2, 3, 4, 0, 2, 1, 4, 0, 2, 7, 3, 4, 2, 7, 6, 4, \
                      2, 5, 1, 4, 2, 5, 6, 4};
  double coord[24], tetCoord[12], centroidTet[3];
  int localIndex[8], localSorted[8];
  for (int i = eptr[e], j = 0; i < eptr[e+1]; ++i, ++j) {
    localIndex[j] = connectivity[i];
    localSorted[j] = connectivity[i];
  }
  qsort(localSorted, 8, sizeof(int), compare);
  for (int i = 0; i < 8; ++i) {
    int indexL = localIndex[i];
    int j = 0;
    while(localSorted[j] != indexL) {
      j = j + 1;
    }
    localIndex[i] = j;
  }
  for (int i = eptr[e], j = 0; i < eptr[e+1]; ++i, ++j) {
    int index0 = ndim*connectivity[i];
    int index1 = ndim*localIndex[j];
    for (int k = 0; k < ndim; ++k) {
      coord[index1+k] = coordinates[index0+k]+displacements[index0+k];
    }
  }
  for (int i = 0; i < 6; ++i) {
    for (int k = 0; k < ndim; ++k) {
      centroidTet[k] = 0.0;
    }
    for (int j = 0; j < 4; ++j) {
      int tetIndex = index[i*4+j]*ndim;
      for (int k = 0; k < ndim; ++k) {
        tetCoord[j*ndim+k] = coord[tetIndex+k];
        centroidTet[k] += tetCoord[j*ndim+k];
      }
    }
    volumeTet = volumeTetrahedron(tetCoord); 
    volume += volumeTet;
    // // Centroid tetrahedron
    // for (int j = 0; j < 4; ++j) {
    //   for (int k = 0; k < ndim; ++k) {
    //     centroidTet[k] += tetCoord[j*ndim+k];
    //   }
    // }
    for (int k = 0; k < ndim; ++k) {
      cm[k] += volumeTet*centroidTet[k]/4.0;
    }
  }
  for (int k = 0; k < ndim; ++k) {
    cm[k] /= volume;
  }
  return volume;
}

double volumeHexahedron(int e) {
  double volume = 0.0, volumeTet;
  // Split hexahedron to 6 tetrahedrons and compute centroid and volume
  // Each set of 4 indices forms a tetrahedron from hexahedron coords
  const int index[24] = {2, 3, 4, 0, 2, 1, 4, 0, 2, 7, 3, 4, 2, 7, 6, 4, \
                      2, 5, 1, 4, 2, 5, 6, 4};
  double coord[24], tetCoord[12];
  int localIndex[8], localSorted[8];
  for (int i = eptr[e], j = 0; i < eptr[e+1]; ++i, ++j) {
    localIndex[j] = connectivity[i];
    localSorted[j] = connectivity[i];
  }
  qsort(localSorted, 8, sizeof(int), compare);
  for (int i = 0; i < 8; ++i) {
    int indexL = localIndex[i];
    int j = 0;
    while(localSorted[j] != indexL) {
      j = j + 1;
    }
    localIndex[i] = j;
  }
  for (int i = eptr[e], j = 0; i < eptr[e+1]; ++i, ++j) {
    int index0 = ndim*connectivity[i];
    int index1 = ndim*localIndex[j];
    for (int k = 0; k < ndim; ++k) {
      coord[index1+k] = coordinates[index0+k];
    }
  }
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 4; ++j) {
      int tetIndex = index[i*4+j]*ndim;
      for (int k = 0; k < ndim; ++k) {
        tetCoord[j*ndim+k] = coord[tetIndex+k];
      }
    }
    volumeTet = volumeTetrahedron(tetCoord); 
    volume += volumeTet;
  }
  return volume;
}
