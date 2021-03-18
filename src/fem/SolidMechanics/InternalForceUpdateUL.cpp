#include "FemTech.h"
#include "blas.h"

void InternalForceUpdateUL(int e, int gp, double *force) {
  const unsigned int start = eptr[e];
  const unsigned int nShapeFunc = eptr[e+1]-start;
  // const unsigned int nShapeFunc = nShapeFunctions[e];
  const unsigned int indexStart = dsptr[e] + gp * GaussPoints[e] * ndim;

  const unsigned int wIndex = gpPtr[e]+gp;
  const double preFactor = gaussWeights[wIndex]*J_Xi;

  for (unsigned int I = 0; I < nShapeFunc; ++I) {
    const unsigned int indexShp = indexStart + I * ndim;
    // Compute dN_I/dx_j = (dN_I/d_Xi)^T F_Xi^-{-1}
    const double dN_Ibydx = dshp[indexShp]*F_XiInverse[0] + 
      dshp[indexShp+1]*F_XiInverse[3] + dshp[indexShp+2]*F_XiInverse[6];
    const double dN_Ibydy = dshp[indexShp]*F_XiInverse[1] + 
      dshp[indexShp+1]*F_XiInverse[4] + dshp[indexShp +2]*F_XiInverse[7];
    const double dN_Ibydz = dshp[indexShp]*F_XiInverse[2] + 
      dshp[indexShp+1]*F_XiInverse[5] + dshp[indexShp +2]*F_XiInverse[8];

    // Compute f_iI = dN_I/dx_j*sigma_{ji}
    // Get location of Cauchy values
    double * sigma_nLocal = &(sigma_n[sigmaptr[e]+6*gp]);

    const double fxI = dN_Ibydx*sigma_nLocal[0]+\
                       dN_Ibydy*sigma_nLocal[5]+dN_Ibydz*sigma_nLocal[4];
    const double fyI = dN_Ibydx*sigma_nLocal[5]+\
                       dN_Ibydy*sigma_nLocal[1]+dN_Ibydz*sigma_nLocal[3];
    const double fzI = dN_Ibydx*sigma_nLocal[4]+\
                       dN_Ibydy*sigma_nLocal[3]+dN_Ibydz*sigma_nLocal[2];
    const unsigned int forceIndex = I*ndim;
    // Add to fint
    force[forceIndex] += preFactor*fxI;
    force[forceIndex+1] += preFactor*fyI;
    force[forceIndex+2] += preFactor*fzI;
  }
	return;
}
