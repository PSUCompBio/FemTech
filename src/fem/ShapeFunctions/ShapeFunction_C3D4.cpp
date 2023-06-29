#include "FemTech.h"

#include <math.h>


void ShapeFunction_C3D4(int e, int gp, double *Chi, double *detJ){
	 //  Purpose : Compute 3-d isoparametric 4 - node tet element shape
	 //            functions and their derivatives w / r x, y, z

	 //  Introduction to Finite Elements in Engineering, 3th Edition
   //  by Chandrupatla and Belegundu

	double chi, eta, iota;
	chi =  Chi[ndim*gp + 0];
	eta =  Chi[ndim*gp + 1];
	iota = Chi[ndim*gp + 2];

	// The shape functions
	int g = GaussPoints[e];
	const int indexShp = gptr[e] + gp * g;
  const int indexDshp = dsptr[e] + gp * g * ndim;

	shp[indexShp + 0] = chi;
	shp[indexShp + 1] = eta;
	shp[indexShp + 2] = iota;
	shp[indexShp + 3] = 1.0 - eta - iota - chi;


	// The first derivatives

	// with respect to chi
	dshp[indexDshp + ndim * 0 + 0] =  1.0;
	dshp[indexDshp + ndim * 1 + 0] =  0.0;
	dshp[indexDshp + ndim * 2 + 0] =  0.0;
	dshp[indexDshp + ndim * 3 + 0] = -1.0;
	// with respect to eta
	dshp[indexDshp + ndim * 0 + 1] =  0.0;
	dshp[indexDshp + ndim * 1 + 1] =  1.0;
	dshp[indexDshp + ndim * 2 + 1] =  0.0;
	dshp[indexDshp + ndim * 3 + 1] = -1.0;
	// with respect to iota
	dshp[indexDshp + ndim * 0 + 2] =  0.0;
	dshp[indexDshp + ndim * 1 + 2] =  0.0;
	dshp[indexDshp + ndim * 2 + 2] =  1.0;
	dshp[indexDshp + ndim * 3 + 2] = -1.0;

	//  Compute jacobian transformation
 	//  Equ 9. 16 3D stress Analysis

	//  x14 y14 z14
	//  x24 y24 z24
	//  x34 y34 z34
	//double xs[ ndim * ndim ];
	double* xs = new double[ndim * ndim];
	int x = 0, y = 1, z = 2;
	int node1 = 0, node2 = 1, node3 = 2, node4 = 3;
	int index=eptr[e];
	// first row
	xs[0]=coordinates[ndim*connectivity[index+node1]+x] - coordinates[ndim*connectivity[index+node4]+x];
	xs[3]=coordinates[ndim*connectivity[index+node1]+y] - coordinates[ndim*connectivity[index+node4]+y];
	xs[6]=coordinates[ndim*connectivity[index+node1]+z] - coordinates[ndim*connectivity[index+node4]+z];
	// second row
	xs[1]=coordinates[ndim*connectivity[index+node2]+x] - coordinates[ndim*connectivity[index+node4]+x];
	xs[4]=coordinates[ndim*connectivity[index+node2]+y] - coordinates[ndim*connectivity[index+node4]+y];
	xs[7]=coordinates[ndim*connectivity[index+node2]+z] - coordinates[ndim*connectivity[index+node4]+z];
	// third row
	xs[2]=coordinates[ndim*connectivity[index+node3]+x] - coordinates[ndim*connectivity[index+node4]+x];
	xs[5]=coordinates[ndim*connectivity[index+node3]+y] - coordinates[ndim*connectivity[index+node4]+y];
	xs[8]=coordinates[ndim*connectivity[index+node3]+z] - coordinates[ndim*connectivity[index+node4]+z];
  double det, J_Inv[9];

  det = inverse3x3Matrix(xs, J_Inv);
  //Chandrupatla Eq. 9.19
  detJ[gp] = fabs(det);

  // Transform derivatives to global co-ordinates
  double c1, c2, c3;
  int baseIndex;
  for (int i = 0; i < nShapeFunctions[e]; ++i) {
    baseIndex = dsptr[e] + gp * g*ndim + ndim * i;
    c1 = dshp[baseIndex]*J_Inv[0]+dshp[baseIndex+1]*J_Inv[3]+dshp[baseIndex+2]*J_Inv[6];
    c2 = dshp[baseIndex]*J_Inv[1]+dshp[baseIndex+1]*J_Inv[4]+dshp[baseIndex+2]*J_Inv[7];
    c3 = dshp[baseIndex]*J_Inv[2]+dshp[baseIndex+1]*J_Inv[5]+dshp[baseIndex+2]*J_Inv[8];

    dshp[baseIndex] = c1;
    dshp[baseIndex+1] = c2;
    dshp[baseIndex+2] = c3;
	delete[]xs;/*Drupal*/
  }
	//for debugging can be removed...
  FILE_LOGMatrix_SINGLE(DEBUGLOGIGNORE, J_Inv, ndim, ndim, "J Inv in e = %d, gp = %d", e, gp);
  return;
}
