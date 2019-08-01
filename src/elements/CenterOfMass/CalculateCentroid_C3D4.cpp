#include "FemTech.h"

void CalculateCentroid_C3D4(int e, double *cm) {
  for (int k = 0; k < ndim; ++k) {
    cm[k] = 0.0;
  }
  for (int i = eptr[e]; i < eptr[e+1]; ++i) {
    int index = ndim*connectivity[i];
    for (int k = 0; k < ndim; ++k) {
      cm[k] += coordinates[index+k]+displacements[index+k];
    }
  }
  for (int k = 0; k < ndim; ++k) {
    cm[k] /= 4.0;
  }
}
