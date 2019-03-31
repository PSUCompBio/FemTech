#include "FemTech.h"

void GetForce_3D() {
  // TODO(Anil) special treatment for first time step
  // Below algorithm works for n > 0
  const int nDOF = nnodes*ndim;
  // Following Belytschko
  // Set force_n to zero
  memset(f_net, 0, nDOF*sizeof(double));
  memset(fi, 0, nDOF*sizeof(double));
  // Loop over elements and Gauss points
	for(int i=0; i<nelements; i++) {
    int nNodes = nShapeFunctions[i];
    // number of shape functions * ndim
    double *fintLocal = (double*)calloc(nNodes*ndim, sizeof(double));
		for(int j=0; j<GaussPoints[i]; j++) {
      // calculate F^n
			CalculateDeformationGradient(i, j);
      // Calculate Determinant of F
			DeterminateF(i, j);
      // Calculate sigma^n
			StressUpdate(i, j);
		  InternalForceUpdate(i, j, fintLocal);
		} //loop on gauss points
    // Move Local internal for to global force
    for (int k = 0; k < nNodes; ++k) {
      int dIndex = connectivity[eptr[i]+k];
      fi[dIndex*ndim+0] += fintLocal[k*ndim+0];
      fi[dIndex*ndim+1] += fintLocal[k*ndim+1];
      fi[dIndex*ndim+2] += fintLocal[k*ndim+2];
    }
    free(fintLocal);
	} // loop on i, nelements
  // Update net force with internal force
  for (int i = 0; i < nnodes*ndim; ++i) {
    f_net[i] -= fi[i];
  }
	return;
}
