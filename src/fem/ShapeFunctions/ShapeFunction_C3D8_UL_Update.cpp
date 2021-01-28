#include "FemTech.h"


void ShapeFunction_C3D8_UL_Update(int e, int gp, double *detJ){
   /* Purpose : Compute 3-d isoparametricderivates of shape function and jacobian at every time step  */

  for (int i = eptr[e]; i < eptr[e + 1]; i++) {
    FILE_LOG_SINGLE(DEBUGLOGIGNORE, "%d %8.4f %8.4f %8.4f", connectivity[i],
      coordinates[ndim*connectivity[i] + 0],
      coordinates[ndim*connectivity[i] + 1],
      coordinates[ndim*connectivity[i] + 2]);
  }

	int g = GaussPoints[e];
	// Compute jacobian transformation
	// If this gets called everytime function is called,
	// I don't like it.
	// The redeclaration needs to find memory and will cut down
	// on speed.
	double xs[ndim*ndim];
  int index = eptr[e];
	for (int j = 0; j < ndim; j++) {
		xs[0+j*ndim] = ((coordinates[ndim*connectivity[index+1] + j]+displacements[ndim*connectivity[index+1] + j]) - (coordinates[ndim*connectivity[index+0] + j]+displacements[ndim*connectivity[index+0] + j])) * dshp[dsptr[e] + gp * g*ndim + ndim * 1 + 0]
		       	 + ((coordinates[ndim*connectivity[index+2] + j]+displacements[ndim*connectivity[index+2] + j]) - (coordinates[ndim*connectivity[index+3] + j]+displacements[ndim*connectivity[index+3] + j])) * dshp[dsptr[e] + gp * g*ndim + ndim * 2 + 0]
			       + ((coordinates[ndim*connectivity[index+5] + j]+displacements[ndim*connectivity[index+5] + j]) - (coordinates[ndim*connectivity[index+4] + j]+displacements[ndim*connectivity[index+4] + j])) * dshp[dsptr[e] + gp * g*ndim + ndim * 5 + 0]
			       + ((coordinates[ndim*connectivity[index+6] + j]+displacements[ndim*connectivity[index+6] + j]) - (coordinates[ndim*connectivity[index+7] + j]+displacements[ndim*connectivity[index+7] + j])) * dshp[dsptr[e] + gp * g*ndim + ndim * 6 + 0];

		xs[1+j*ndim] = ((coordinates[ndim*connectivity[index+2] + j]+displacements[ndim*connectivity[index+2] + j]) - (coordinates[ndim*connectivity[index+1] + j]+displacements[ndim*connectivity[index+1] + j])) * dshp[dsptr[e] + gp * g*ndim + ndim * 2 + 1]
			       + ((coordinates[ndim*connectivity[index+3] + j]+displacements[ndim*connectivity[index+3] + j]) - (coordinates[ndim*connectivity[index+0] + j]+displacements[ndim*connectivity[index+0] + j])) * dshp[dsptr[e] + gp * g*ndim + ndim * 3 + 1]
			       + ((coordinates[ndim*connectivity[index+6] + j]+displacements[ndim*connectivity[index+6] + j]) - (coordinates[ndim*connectivity[index+5] + j]+displacements[ndim*connectivity[index+5] + j])) * dshp[dsptr[e] + gp * g*ndim + ndim * 6 + 1]
			       + ((coordinates[ndim*connectivity[index+7] + j]+displacements[ndim*connectivity[index+7] + j]) - (coordinates[ndim*connectivity[index+4] + j]+displacements[ndim*connectivity[index+4] + j])) * dshp[dsptr[e] + gp * g*ndim + ndim * 7 + 1];

		xs[2+j*ndim] = ((coordinates[ndim*connectivity[index+4] + j]+displacements[ndim*connectivity[index+4] + j]) - (coordinates[ndim*connectivity[index+0] + j]+displacements[ndim*connectivity[index+0] + j])) * dshp[dsptr[e] + gp * g*ndim + ndim * 4 + 2]
		         + ((coordinates[ndim*connectivity[index+5] + j]+displacements[ndim*connectivity[index+5] + j]) - (coordinates[ndim*connectivity[index+1] + j]+displacements[ndim*connectivity[index+1] + j])) * dshp[dsptr[e] + gp * g*ndim + ndim * 5 + 2]
			       + ((coordinates[ndim*connectivity[index+6] + j]+displacements[ndim*connectivity[index+6] + j]) - (coordinates[ndim*connectivity[index+2] + j]+displacements[ndim*connectivity[index+2] + j])) * dshp[dsptr[e] + gp * g*ndim + ndim * 6 + 2]
			       + ((coordinates[ndim*connectivity[index+7] + j]+displacements[ndim*connectivity[index+7] + j]) - (coordinates[ndim*connectivity[index+3] + j]+displacements[ndim*connectivity[index+3] + j])) * dshp[dsptr[e] + gp * g*ndim + ndim * 7 + 2];
	}

  double det, J_Inv[9];

  inverse3x3Matrix(xs, J_Inv, &det);
  detJ[gp] = det;

  // Transform derivatives to global co-ordinates
  double c1, c2, c3;
  int baseIndex;

//	FILE_LOG_SINGLE(DEBUGLOGIGNORE, "Shape Function C3D8\nDeterminant of Jacobian : "
//      "%12.6e\nDerivatives eid : %d, gpid : %d, chi : %12.6f, eta : %12.6f, "
//      "iota : %12.6f", det, e, gp, chi, eta, iota);

  for (int i = 0; i < nShapeFunctions[e]; ++i) {
    baseIndex = dsptr[e] + gp * g*ndim + ndim * i;
    c1 = dshp[baseIndex]*J_Inv[0]+dshp[baseIndex+1]*J_Inv[3]+dshp[baseIndex+2]*J_Inv[6];
    c2 = dshp[baseIndex]*J_Inv[1]+dshp[baseIndex+1]*J_Inv[4]+dshp[baseIndex+2]*J_Inv[7];
    c3 = dshp[baseIndex]*J_Inv[2]+dshp[baseIndex+1]*J_Inv[5]+dshp[baseIndex+2]*J_Inv[8];

	  FILE_LOG_SINGLE(DEBUGLOGIGNORE, "Shape fn        : %d, %12.6f, %12.6f, %12.6f", i, \
        dshp[baseIndex], dshp[baseIndex+1], dshp[baseIndex+2]);

    dshp[baseIndex] = c1;
    dshp[baseIndex+1] = c2;
    dshp[baseIndex+2] = c3;

    FILE_LOG_SINGLE(DEBUGLOGIGNORE, "Shape fn Global : %d, %12.6f, %12.6f, %12.6f", \
        i, dshp[baseIndex], dshp[baseIndex+1], dshp[baseIndex+2]);
  }

  FILE_LOGMatrix_SINGLE(DEBUGLOGIGNORE, xs, ndim, ndim, "---- XS Matrix ----")
  FILE_LOGMatrix_SINGLE(DEBUGLOGIGNORE, J_Inv, ndim, ndim, "---- JInv Matrix ----")
  return;
}
