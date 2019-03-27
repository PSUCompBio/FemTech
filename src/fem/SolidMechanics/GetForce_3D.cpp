#include "FemTech.h"
#include "blas.h"

void GetForce_3D() {
  // TODO(Anil) special treatment for first time step
  // Below algorithm works for n > 0
  const int nDOF = nnodes*ndim;
  // Following Belytschko
  // Set force_n to zero
  memset(f_net, 0, nDOF*sizeof(double));
  int cSize = 6;
  // Variables to store intermediate outputs
  double *Bdn = (double*)calloc(cSize, sizeof(double));
  double *CBdn = (double*)calloc(cSize, sizeof(double));
  // Loop over elements and Gauss points 
	for(int i=0; i<nelements; i++) {
    int nNodes = nShapeFunctions[i];
    // number of shape functions * ndim
    int bColSize = nNodes*ndim;
    int Bsize = bColSize*cSize;
    double *B = (double*)calloc(Bsize, sizeof(double));
    double *localDisplacement = (double*)calloc(nNodes*ndim, sizeof(double));
    double *fintLocal = (double*)calloc(nNodes*ndim, sizeof(double));
    // Copy local displacement to array
    for (int k = 0; k < nNodes; ++k) {
      int dIndex = connectivity[eptr[i]+k];
      localDisplacement[k*ndim+0] = displacements[dIndex*ndim+0];
      localDisplacement[k*ndim+1] = displacements[dIndex*ndim+1];
      localDisplacement[k*ndim+2] = displacements[dIndex*ndim+2];
    }
		for(int j=0; j<GaussPoints[i]; j++) {
      // calculate F^n
			CalculateDeformationGradient(i, j);
      // Calculate Determinant of F
			DeterminateF(i, j);
      // Calculate B matrix for each shape function
      for (int k = 0; k < nShapeFunctions[i]; ++k) {
        StressDisplacementMatrix(i, j, k, &(B[6*ndim*k]));
      }
      // Compute B*d^n
      dgemv_(chn, &cSize, &bColSize, &one, B, &cSize, localDisplacement, \
          &oneI, &zero, Bdn, &oneI);
      // Compute sigma^n = C*B*d^n
      dgemv_(chn, &cSize, &cSize, &one, C, &cSize, Bdn, \
          &oneI, &zero, CBdn, &oneI);
      // Compute B^T*sigma^n
      dgemv_(chy, &cSize, &bColSize, &one, B, &cSize, CBdn, \
          &oneI, &zero, fintLocal, &oneI);
      // Add to fint
      int wIndex = gpPtr[i]+j;
      const double preFactor = gaussWeights[wIndex]*detF[wIndex];
      for (int k = 0; k < bColSize; ++k) {
        fintLocal[k] += preFactor*fintLocal[k];
      }
			// InverseF(i,j);
			// StressUpdate(i,j);
			// InternalForceUpdate(i,j);
		} //loop on gauss points
    for (int k = 0; k < nNodes; ++k) {
      int dIndex = connectivity[eptr[i]+k];
      f_net[dIndex*ndim+0] -= fintLocal[k*ndim+0];
      f_net[dIndex*ndim+1] -= fintLocal[k*ndim+1];
      f_net[dIndex*ndim+2] -= fintLocal[k*ndim+2];
    }
    free(B);
    free(localDisplacement);
    free(fintLocal);
	} // loop on i, nelements
  free(Bdn);
  free(CBdn);
	return;
}
