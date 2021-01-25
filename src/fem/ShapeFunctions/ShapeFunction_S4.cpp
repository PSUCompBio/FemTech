#include "FemTech.h"


void ShapeFunction_S4(int e, int gp, double *Chi, double *detJ){
   /* Purpose : Compute 3-d isoparametric 8 - node e shape
	            functions and their derivatives with respect to
	            chi, eta, iota

				 NOTE: Might need to add dirivatives with respect to
					   x, y, z
  */

	double chi, eta;
	chi =  Chi[ndim*gp + 0];
	eta =  Chi[ndim*gp + 1];
	ndim = 3;
	int ndimshell = 2;

	// The shape functions
	int g = GaussPoints[e];
  const int indexShp = gptr[e] + gp * g;
  const int indexDshp = dsptr[e] + gp * g * ndimshell;


	shp[indexShp + 0] = ((1 - chi)*(1 - eta)) / 4;
	shp[indexShp + 1] = ((1 + chi)*(1 - eta)) / 4;
	shp[indexShp + 2] = ((1 + chi)*(1 + eta)) / 4;
	shp[indexShp + 3] = ((1 - chi)*(1 + eta)) / 4;


	// The first derivatives

  // with respect to chi
	dshp[indexDshp + ndimshell * 0 + 0] = -(1 - eta) / 4;
	dshp[indexDshp + ndimshell * 1 + 0] =  (1 - eta) / 4;
	dshp[indexDshp + ndimshell * 2 + 0] =  (1 + eta) / 4;
	dshp[indexDshp + ndimshell * 3 + 0] = -(1 + eta) / 4;

	// with respect to eta
	dshp[indexDshp + ndimshell * 0 + 1] = -(1 - chi) / 4;
	dshp[indexDshp + ndimshell * 1 + 1] = -(1 + chi) / 4;
	dshp[indexDshp + ndimshell * 2 + 1] =  (1 + chi) / 4;
	dshp[indexDshp + ndimshell * 3 + 1] =  (1 - chi) / 4;


  for (int i = eptr[e]; i < eptr[e + 1]; i++) {
    FILE_LOG_SINGLE(DEBUGLOGIGNORE, "%d %8.4f %8.4f", connectivity[i],
      coordinates[ndim*connectivity[i] + 0],
      coordinates[ndim*connectivity[i] + 1]);
  }
	// Compute jacobian transformation
	// If this gets called everytime function is called,
	// I don't like it.
	// The redeclaration needs to find memory and will cut down
	// on speed.
	double xs[ndimshell*ndimshell];
  int index = eptr[e];
	for (int j = 0; j < ndimshell; j++) {
		xs[0+j*ndimshell] = (hatcoordinates[ndim*connectivity[index+1] + j] - hatcoordinates[ndim*connectivity[index+0] + j]) * dshp[dsptr[e] + gp * g*ndim + ndim * 1 + 0]
		       	 + (hatcoordinates[ndim*connectivity[index+2] + j] - hatcoordinates[ndim*connectivity[index+3] + j]) * dshp[dsptr[e] + gp * g*ndim + ndim * 2 + 0];

		xs[1+j*ndimshell] = (hatcoordinates[ndim*connectivity[index+2] + j] - hatcoordinates[ndim*connectivity[index+1] + j]) * dshp[dsptr[e] + gp * g*ndimshell + ndimshell * 2 + 1]
			       + (hatcoordinates[ndim*connectivity[index+3] + j] - hatcoordinates[ndim*connectivity[index+0] + j]) * dshp[dsptr[e] + gp * g*ndimshell + ndimshell * 3 + 1];	
	}

	double det, J_Inv[4];
  det = xs[0]*xs[3]-xs[1]*xs[2];
  J_Inv[0] = xs[3]/det;
	J_Inv[1] = -xs[1]/det;
	J_Inv[2] = -xs[2]/det;
	J_Inv[3] = xs[0]/det;
  detJ[gp] = fabs(det);

  // Transform derivatives to global co-ordinates
  double c1, c2;
  int baseIndex;

	FILE_LOG_SINGLE(DEBUGLOGIGNORE, "Shape Function S4\nDeterminant of Jacobian : "
      "%12.6e\nDerivatives eid : %d, gpid : %d, chi : %12.6f, eta : %12.6f, "
      , det, e, gp, chi, eta);

  for (int i = 0; i < nShapeFunctions[e]; ++i) {
    baseIndex = dsptr[e] + gp * g*ndim + ndim * i;
    c1 = dshp[baseIndex]*J_Inv[0]+dshp[baseIndex+1]*J_Inv[2];
    c2 = dshp[baseIndex]*J_Inv[1]+dshp[baseIndex+1]*J_Inv[3];

	  FILE_LOG_SINGLE(DEBUGLOGIGNORE, "Shape fn        : %d, %12.6f, %12.6f", i, \
        dshp[baseIndex], dshp[baseIndex+1]);

    dshp[baseIndex] = c1;
    dshp[baseIndex+1] = c2;

    FILE_LOG_SINGLE(DEBUGLOGIGNORE, "Shape fn Global : %d, %12.6f, %12.6f", \
        i, dshp[baseIndex], dshp[baseIndex+1]);
  }

  FILE_LOGMatrix_SINGLE(DEBUGLOGIGNORE, xs, ndim, ndim, "---- XS Matrix ----")
  FILE_LOGMatrix_SINGLE(DEBUGLOGIGNORE, J_Inv, ndim, ndim, "---- JInv Matrix ----")
  return;
}
