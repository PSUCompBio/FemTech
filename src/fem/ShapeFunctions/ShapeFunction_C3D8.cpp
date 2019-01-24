#include "FemTech.h"


void ShapeFunction_C3D8(int e, int gp, double *Chi){
   /*
	   subroutine shp3d(ss, xsj, shp, xl, ndm)

	   Purpose : Compute 3-d isoparametric 8 - node e shape
	             functions and their derivatives with respect to 
				 chi, eta, iota 

				 NOTE: Might need to add dirivatives with respect to
					   x, y, z	

	   Inputs :
	   Chi- Natural coordinates of point
	   xl(ndm, *) - Nodal coordinates for e
	   ndim - Spatial dimension of mesh

	   Outputs :
	   xsj - Jacobian determinant at point
	   shp(4, *) - Shape functions and derivatives at point
	   shp(1, i) = dN_i / dx
	   shp(2, i) = dN_i / dy
	   shp(3, i) = dN_i / dz
	   shp(4, i) = N_i
	  
	   */
	int debug = 1;

	double chi, eta, iota;
	chi =  Chi[ndim*gp + 0];
	eta =  Chi[ndim*gp + 1];
	iota = Chi[ndim*gp + 2];

//	printf("chi = %f\neta = %f\niota = %f | %f\n", chi, eta, iota,
//		-((chi - 1)*(eta - 1)*(iota - 1)) / 8);

	// The shape functions
	int g = GaussPoints[e];
	shp[gptr[e] + gp * g + 0] = -((chi - 1)*(eta - 1)*(iota - 1)) / 8;
	shp[gptr[e] + gp * g + 1] =  ((chi + 1)*(eta - 1)*(iota - 1)) / 8;
	shp[gptr[e] + gp * g + 2] =  ((chi - 1)*(eta + 1)*(iota - 1)) / 8;
	shp[gptr[e] + gp * g + 3] = -((chi + 1)*(eta + 1)*(iota - 1)) / 8;
	shp[gptr[e] + gp * g + 4] =  ((chi - 1)*(eta - 1)*(iota + 1)) / 8;
	shp[gptr[e] + gp * g + 5] = -((chi + 1)*(eta - 1)*(iota + 1)) / 8;
	shp[gptr[e] + gp * g + 6] = -((chi - 1)*(eta + 1)*(iota + 1)) / 8;
	shp[gptr[e] + gp * g + 7] =  ((chi + 1)*(eta + 1)*(iota + 1)) / 8;

	if (debug) {
		printf("e.%d: ", e);
		//for (int i = 0; i < g; i++) {
		//	printf("%f ", shp[gptr[e] + gp * g + i]);
		//}
		

		printf("%f ", -((chi - 1)*(eta - 1)*(iota - 1)) / 8);
		printf("%f ", ((chi + 1)*(eta - 1)*(iota - 1)) / 8);
		printf("%f ", ((chi - 1)*(eta + 1)*(iota - 1)) / 8);
		printf("%f ", -((chi + 1)*(eta + 1)*(iota - 1)) / 8);
		printf("%f ", ((chi - 1)*(eta - 1)*(iota + 1)) / 8);
		printf("%f ", -((chi + 1)*(eta - 1)*(iota + 1)) / 8);
		printf("%f ", -((chi - 1)*(eta + 1)*(iota + 1)) / 8);
		printf("%f ", ((chi + 1)*(eta + 1)*(iota + 1)) / 8);
		printf("\n");
	}

	// The first derivatives

    // with respect to chi
	dshp[dsptr[e] + gp*g*ndim + ndim*0 + 0] = -((eta - 1)*(iota - 1)) / 8;
	dshp[dsptr[e] + gp*g*ndim + ndim*1 + 0] =  ((eta - 1)*(iota - 1)) / 8;
	dshp[dsptr[e] + gp*g*ndim + ndim*2 + 0] =  ((eta + 1)*(iota - 1)) / 8;
	dshp[dsptr[e] + gp*g*ndim + ndim*3 + 0] = -((eta + 1)*(iota - 1)) / 8;
	dshp[dsptr[e] + gp*g*ndim + ndim*4 + 0] =  ((eta - 1)*(iota + 1)) / 8;
	dshp[dsptr[e] + gp*g*ndim + ndim*5 + 0] = -((eta - 1)*(iota + 1)) / 8;
	dshp[dsptr[e] + gp*g*ndim + ndim*6 + 0] = -((eta + 1)*(iota + 1)) / 8;
	dshp[dsptr[e] + gp*g*ndim + ndim*7 + 0] =  ((eta + 1)*(iota + 1)) / 8;

	// with respect to eta
	dshp[dsptr[e] + gp * g*ndim + ndim * 0 + 1] = -((chi - 1)*(iota - 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 1 + 1] =  ((chi + 1)*(iota - 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 2 + 1] =  ((chi - 1)*(iota - 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 3 + 1] = -((chi + 1)*(iota - 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 4 + 1] =  ((chi - 1)*(iota + 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 5 + 1] = -((chi + 1)*(iota + 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 6 + 1] = -((chi - 1)*(iota + 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 7 + 1] =  ((chi + 1)*(iota + 1)) / 8;

	// with respect to  iota
	dshp[dsptr[e] + gp * g*ndim + ndim * 0 + 2] = -((chi - 1)*(eta - 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 1 + 2] =  ((chi + 1)*(eta - 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 2 + 2] =  ((chi - 1)*(eta + 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 3 + 2] = -((chi + 1)*(eta + 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 4 + 2] =  ((chi - 1)*(eta - 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 5 + 2] = -((chi + 1)*(eta - 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 6 + 2] = -((chi - 1)*(eta + 1)) / 8;
	dshp[dsptr[e] + gp * g*ndim + ndim * 7 + 2] =  ((chi + 1)*(eta + 1)) / 8;
	  

   return;
}
