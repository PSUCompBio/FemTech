#include "FemTech.h"


void ShapeFunction_C3D8(int e, int gp, double *Chi){
   /* Purpose : Compute 3-d isoparametric 8 - node e shape
	            functions and their derivatives with respect to 
	            chi, eta, iota 

				 NOTE: Might need to add dirivatives with respect to
					   x, y, z	
  */
	int debug = 0;

	double chi, eta, iota;
	chi =  Chi[ndim*gp + 0];
	eta =  Chi[ndim*gp + 1];
	iota = Chi[ndim*gp + 2];

	// The shape functions
	int g = GaussPoints[e];

	shp[gptr[e] + gp * g + 0] = ((1 - chi)*(1 - eta)*(1 - iota)) / 8;
	shp[gptr[e] + gp * g + 1] = ((1 + chi)*(1 - eta)*(1 - iota)) / 8;
	shp[gptr[e] + gp * g + 2] = ((1 + chi)*(1 + eta)*(1 - iota)) / 8;
	shp[gptr[e] + gp * g + 3] = ((1 - chi)*(1 + eta)*(1 - iota)) / 8;
	shp[gptr[e] + gp * g + 4] = ((1 - chi)*(1 - eta)*(1 + iota)) / 8;
	shp[gptr[e] + gp * g + 5] = ((1 + chi)*(1 - eta)*(1 + iota)) / 8;
	shp[gptr[e] + gp * g + 6] = ((1 + chi)*(1 + eta)*(1 + iota)) / 8;
	shp[gptr[e] + gp * g + 7] = ((1 - chi)*(1 + eta)*(1 + iota)) / 8;


	// The first derivatives

    // with respect to chi
	dshp[dsptr[e] + gp * g*ndim + ndim * 0 + 0] = -((eta - 1)*(iota - 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 1 + 0] =  ((eta - 1)*(iota - 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 2 + 0] = -((eta + 1)*(iota - 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 3 + 0] =  ((eta + 1)*(iota - 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 4 + 0] =  ((eta - 1)*(iota + 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 5 + 0] = -((eta - 1)*(iota + 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 6 + 0] =  ((eta + 1)*(iota + 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 7 + 0] = -((eta + 1)*(iota + 1)) / 8;

	// with respect to eta
	dshp[dsptr[e] + gp * g*ndim + ndim * 0 + 1] = -((chi - 1)*(iota - 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 1 + 1] =  ((chi + 1)*(iota - 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 2 + 1] = -((chi + 1)*(iota - 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 3 + 1] =  ((chi - 1)*(iota - 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 4 + 1] =  ((chi - 1)*(iota + 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 5 + 1] = -((chi + 1)*(iota + 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 6 + 1] =  ((chi + 1)*(iota + 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 7 + 1] = -((chi - 1)*(iota + 1)) / 8;

	// with respect to  iota
	dshp[dsptr[e] + gp * g*ndim + ndim * 0 + 2] = -((chi - 1)*(eta - 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 1 + 2] =  ((chi + 1)*(eta - 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 2 + 2] = -((chi + 1)*(eta + 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 3 + 2] =  ((chi - 1)*(eta + 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 4 + 2] =  ((chi - 1)*(eta - 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 5 + 2] = -((chi + 1)*(eta - 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 6 + 2] =  ((chi + 1)*(eta + 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 7 + 2] = -((chi - 1)*(eta + 1)) / 8;


	//for debugging can be removed...
	if (debug == 1 && 1 == 0) {
		//printf("e.%d ", e);
		for (int i = eptr[e]; i < eptr[e + 1]; i++) {
			//printf("%d ", i);
			//printf("%d ", connectivity[i]);
			printf("%d %8.4f %8.4f %8.4f \n", connectivity[i],
				coordinates[ndim*connectivity[i] + 0],
				coordinates[ndim*connectivity[i] + 1],
				coordinates[ndim*connectivity[i] + 2]);
		}
		printf("\n");
	}


	// Compute jacobian transformation

	// R. Kraft does not like using multidimensional arrays. 
	// this wss done just to keep making prgress
	// ALSO, if this gets called everytime function is called, 
	// I don't like it.
	// The redeclaration needs to find memory and will cut down
	// on speed. 
	double xs[3][3]; 
	for (int j = 0; j < 3; j++) {
		xs[0][j] = (coordinates[ndim*connectivity[1] + j] - coordinates[ndim*connectivity[0] + j]) * dshp[dsptr[e] + gp * g*ndim + ndim * 1 + 0]
		       	 + (coordinates[ndim*connectivity[2] + j] - coordinates[ndim*connectivity[3] + j]) * dshp[dsptr[e] + gp * g*ndim + ndim * 2 + 0]
			     + (coordinates[ndim*connectivity[5] + j] - coordinates[ndim*connectivity[4] + j]) * dshp[dsptr[e] + gp * g*ndim + ndim * 5 + 0]
			     + (coordinates[ndim*connectivity[6] + j] - coordinates[ndim*connectivity[7] + j]) * dshp[dsptr[e] + gp * g*ndim + ndim * 6 + 0];

		xs[1][j] = (coordinates[ndim*connectivity[2] + j] - coordinates[ndim*connectivity[1] + j]) * dshp[dsptr[e] + gp * g*ndim + ndim * 2 + 1]
			     + (coordinates[ndim*connectivity[3] + j] - coordinates[ndim*connectivity[0] + j]) * dshp[dsptr[e] + gp * g*ndim + ndim * 3 + 1]
			     + (coordinates[ndim*connectivity[6] + j] - coordinates[ndim*connectivity[5] + j]) * dshp[dsptr[e] + gp * g*ndim + ndim * 6 + 1]
			     + (coordinates[ndim*connectivity[7] + j] - coordinates[ndim*connectivity[4] + j]) * dshp[dsptr[e] + gp * g*ndim + ndim * 7 + 1];

		xs[2][j] = (coordinates[ndim*connectivity[4] + j] - coordinates[ndim*connectivity[0] + j]) * dshp[dsptr[e] + gp * g*ndim + ndim * 4 + 2]
		         + (coordinates[ndim*connectivity[5] + j] - coordinates[ndim*connectivity[1] + j]) * dshp[dsptr[e] + gp * g*ndim + ndim * 5 + 2]
			     + (coordinates[ndim*connectivity[6] + j] - coordinates[ndim*connectivity[2] + j]) * dshp[dsptr[e] + gp * g*ndim + ndim * 6 + 2]
			     + (coordinates[ndim*connectivity[7] + j] - coordinates[ndim*connectivity[3] + j]) * dshp[dsptr[e] + gp * g*ndim + ndim * 7 + 2];
	}

	//for debugging can be removed...
	if (debug == 1 && 1 == 0) {
		printf("%8.4e %8.4e %8.4e\n", xs[0][0], xs[0][1], xs[0][2]);
		printf("%8.4e %8.4e %8.4e\n", xs[1][0], xs[1][1], xs[1][2]);
		printf("%8.4e %8.4e %8.4e\n", xs[2][0], xs[2][1], xs[2][2]);
		printf("\n");
	}
   return;
}
