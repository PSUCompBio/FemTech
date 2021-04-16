#include "FemTech.h"
#include "blas.h"

void InternalForceUpdateUL(int e, int gp, double *force) {
  const unsigned int start = eptr[e];
  const unsigned int nShapeFunc = eptr[e+1]-start;
  // const unsigned int nShapeFunc = nShapeFunctions[e];
  const unsigned int indexStart = dsptr[e] + gp * GaussPoints[e] * ndim;

  const unsigned int wIndex = gpPtr[e]+gp;
  const double preFactor = gaussWeights[wIndex]*J_Xi;
//  printf("%.10f %.10f\n", preFactor, Time);
  double p, p1, p2;
  double rho = properties[MAXMATPARAMS * pid[e] + 0];
  double mu = properties[MAXMATPARAMS * pid[e] + 1];
  double lambda = properties[MAXMATPARAMS * pid[e] + 2];
  double ce = sqrt((lambda+2*mu)/rho);
  double le = CalculateCharacteristicLength(e);
  double strate = 0.0;
  double b1=0.06;
  double b2=1.2;
  for(int i = eptr[e], I = 0; i < eptr[e+1]; ++i, ++I)
  {
    const unsigned int indexShp = indexStart + I * ndim;
    strate += velocities_half[ndim*connectivity[i]+0]*(dshp[indexShp]*F_XiInverse[0] + dshp[indexShp+1]*F_XiInverse[3] + dshp[indexShp+2]*F_XiInverse[6]) + velocities_half[ndim*connectivity[i]+1]*(dshp[indexShp]*F_XiInverse[1] + dshp[indexShp+1]*F_XiInverse[4] + dshp[indexShp +2]*F_XiInverse[7]) + velocities_half[ndim*connectivity[i]+2]*(dshp[indexShp]*F_XiInverse[2] + dshp[indexShp+1]*F_XiInverse[5] + dshp[indexShp +2]*F_XiInverse[8]);
    }
    printf("%.10f %.10f %.10f %.10f\n", velocities_half[ndim*3+0], velocities_half[ndim*3+1], velocities_half[ndim*3+2], Time);
  p1 = b1*rho*ce*le*strate;
  p2 = rho*pow((b2*le),2)*strate*std::min(0.0, strate);
  p = 0.0;
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
    const double fxI = dN_Ibydx*(sigma_n[0]+p)+dN_Ibydy*(sigma_n[5]+p)+dN_Ibydz*(sigma_n[4]+p);
    const double fyI = dN_Ibydx*(sigma_n[5]+p)+dN_Ibydy*(sigma_n[1]+p)+dN_Ibydz*(sigma_n[3]+p);
    const double fzI = dN_Ibydx*(sigma_n[4]+p)+dN_Ibydy*(sigma_n[3]+p)+dN_Ibydz*(sigma_n[2]+p);
    const unsigned int forceIndex = I*ndim;
    // Add to fint
    force[forceIndex] += preFactor*fxI;
    force[forceIndex+1] += preFactor*fyI;
    force[forceIndex+2] += preFactor*fzI;
  }
  return;
}
