#include "FemTech.h"


void ShapeFunction_C3D8(int element, int intpt, int nGaussPoints, double *Chi){
   /*
	   subroutine shp3d(ss, xsj, shp, xl, ndm)

	   Purpose : Compute 3-d isoparametric 8 - node element shape
	             functions and their derivatives with respect to 
				 chi, eta, iota 

				 NOTE: Might need to add dirivatives with respect to
					   x, y, z	

	   Inputs :
	   Chi- Natural coordinates of point
	   xl(ndm, *) - Nodal coordinates for element
	   ndim - Spatial dimension of mesh

	   Outputs :
	   xsj - Jacobian determinant at point
	   shp(4, *) - Shape functions and derivatives at point
	   shp(1, i) = dN_i / dx
	   shp(2, i) = dN_i / dy
	   shp(3, i) = dN_i / dz
	   shp(4, i) = N_i
	  
	   */

	double chi, eta, iota;
	chi = Chi[ndim*intpt + 0];
	eta = Chi[ndim*intpt + 1];
	iota = Chi[ndim*intpt + 2];

	// The shape functions
	shp[gptr[element] + 0] = -((chi - 1)*(eta - 1)*(iota - 1)) / 8;
	shp[gptr[element] + 1] =  ((chi + 1)*(eta - 1)*(iota - 1)) / 8;
	shp[gptr[element] + 2] =  ((chi - 1)*(eta + 1)*(iota - 1)) / 8;
	shp[gptr[element] + 3] = -((chi + 1)*(eta + 1)*(iota - 1)) / 8;
	shp[gptr[element] + 4] =  ((chi - 1)*(eta - 1)*(iota + 1)) / 8;
	shp[gptr[element] + 5] = -((chi + 1)*(eta - 1)*(iota + 1)) / 8;
	shp[gptr[element] + 6] = -((chi - 1)*(eta + 1)*(iota + 1)) / 8;
	shp[gptr[element] + 7] =  ((chi + 1)*(eta + 1)*(iota + 1)) / 8;

	// The first derivatives

    // with respect to chi
	dshp[dsptr[element] + 0 * ndim + 0] = -((eta - 1)*(iota - 1)) / 8;
	dshp[dsptr[element] + 1 * ndim + 0] =  ((eta - 1)*(iota - 1)) / 8;
	dshp[dsptr[element] + 2 * ndim + 0] =  ((eta + 1)*(iota - 1)) / 8;
	dshp[dsptr[element] + 3 * ndim + 0] = -((eta + 1)*(iota - 1)) / 8;
	dshp[dsptr[element] + 4 * ndim + 0] =  ((eta - 1)*(iota + 1)) / 8;
	dshp[dsptr[element] + 5 * ndim + 0] = -((eta - 1)*(iota + 1)) / 8;
	dshp[dsptr[element] + 6 * ndim + 0] = -((eta + 1)*(iota + 1)) / 8;
	dshp[dsptr[element] + 7 * ndim + 0] =  ((eta + 1)*(iota + 1)) / 8;

	// with respect to eta
	dshp[dsptr[element] + 0 * ndim + 1] = -((chi - 1)*(iota - 1)) / 8;
	dshp[dsptr[element] + 1 * ndim + 1] =  ((chi + 1)*(iota - 1)) / 8;
	dshp[dsptr[element] + 2 * ndim + 1] =  ((chi - 1)*(iota - 1)) / 8;
	dshp[dsptr[element] + 3 * ndim + 1] = -((chi + 1)*(iota - 1)) / 8;
	dshp[dsptr[element] + 4 * ndim + 1] =  ((chi - 1)*(iota + 1)) / 8;
	dshp[dsptr[element] + 5 * ndim + 1] = -((chi + 1)*(iota + 1)) / 8;
	dshp[dsptr[element] + 6 * ndim + 1] = -((chi - 1)*(iota + 1)) / 8;
	dshp[dsptr[element] + 7 * ndim + 1] =  ((chi + 1)*(iota + 1)) / 8;

	// with respect to  iota
	dshp[dsptr[element] + 0 * ndim + 2] = -((chi - 1)*(eta - 1)) / 8;
	dshp[dsptr[element] + 1 * ndim + 2] =  ((chi + 1)*(eta - 1)) / 8;
	dshp[dsptr[element] + 2 * ndim + 2] =  ((chi - 1)*(eta + 1)) / 8;
	dshp[dsptr[element] + 3 * ndim + 2] = -((chi + 1)*(eta + 1)) / 8;
	dshp[dsptr[element] + 4 * ndim + 2] =  ((chi - 1)*(eta - 1)) / 8;
	dshp[dsptr[element] + 5 * ndim + 2] = -((chi + 1)*(eta - 1)) / 8;
	dshp[dsptr[element] + 6 * ndim + 2] = -((chi - 1)*(eta + 1)) / 8;
	dshp[dsptr[element] + 7 * ndim + 2] =  ((chi + 1)*(eta + 1)) / 8;
	  
/*

	   //integer   ndm, i, j, k
   //int   ndm, i, j, k;
	   //real * 8    rxsj, xsj, ap1, am1, ap2, am2, ap3, am3, c1, c2, c3
   //double    rxsj, xsj, ap1, am1, ap2, am2, ap3, am3, c1, c2, c3;

   //real * 8    ss(3), shp(4, 8), xl(ndm, 8), xs(3, 3), ad(3, 3)
	   //double    ss[3-1], shp[4-1, 8-1], xl[ndm-1, 8-1], xs[3-1, 3-1], ad[3-1, 3-1];

   /*
   double shp[4 - 1][8 - 1], xs[3 - 1][3 - 1];

	   //save

	   //!Compute shape functions and their natural coord.derivatives



	   double ap1 = 1.0 + ss[1 - 1];
	   double am1 = 1.0 - ss[1 - 1];
	   double ap2 = 1.0 + ss[2 - 1];
	   double am2 = 1.0 - ss[2 - 1];
	   double ap3 = 1.0 + ss[3 - 1];
	   double am3 = 1.0 - ss[3 - 1];

	   //!Compute for (-, -) values

	   double c1 = 0.125*am1*am2;
	   double c2 = 0.125*am2*am3;
	   double c3 = 0.125*am1*am3;
	   shp[1 - 1][1 - 1] = -c2;
	   shp[1 - 1][2 - 1] = c2;
	   shp[2 - 1][1 - 1] = -c3;
	   shp[2 - 1][4 - 1] = c3;
	   shp[3 - 1][1 - 1] = -c1;
	   shp[3 - 1][5 - 1] = c1;
	   shp[4 - 1][1 - 1] = c1 * am3;
	   shp[4 - 1][5 - 1] = c1 * ap3;

	   //!Compute for (+, +) values

	   c1 = 0.125*ap1*ap2;
	   c2 = 0.125*ap2*ap3;
	   c3 = 0.125*ap1*ap3;
	   shp[1-1][8 - 1] = -c2;
	   shp[1 - 1][7 - 1] = c2;
	   shp[2 - 1][6 - 1] = -c3;
	   shp[2 - 1][7 - 1] = c3;
	   shp[3 - 1][3 - 1] = -c1;
	   shp[3 - 1][7 - 1] = c1;
	   shp[4 - 1][3 - 1] = c1 * am3;
	   shp[4 - 1][7 - 1] = c1 * ap3;

	   //!Compute for (-, +) values

	   c1 = 0.125*am1*ap2;
	   c2 = 0.125*am2*ap3;
	   c3 = 0.125*am1*ap3;
	   shp[1 - 1][5 - 1] = -c2;
	   shp[1 - 1][6 - 1] = c2;
	   shp[2 - 1][5 - 1] = -c3;
	   shp[2 - 1][8 - 1] = c3;
	   shp[3 - 1][4 - 1] = -c1;
	   shp[3 - 1][8 - 1] = c1;
	   shp[4 - 1][4 - 1] = c1 * am3;
	   shp[4 - 1][8 - 1] = c1 * ap3;

	   //!Compute for (+, -) values

	   c1 = 0.125*ap1*am2;
	   c2 = 0.125*ap2*am3;
	   c3 = 0.125*ap1*am3;
	   shp[1 - 1][4 - 1] = -c2;
	   shp[1 - 1][3 - 1] = c2;
	   shp[2 - 1][2 - 1] = -c3;
	   shp[2 - 1][3 - 1] = c3;
	   shp[3 - 1][2 - 1] = -c1;
	   shp[3 - 1][6 - 1] = c1;
	   shp[4 - 1][2 - 1] = c1 * am3;
	   shp[4 - 1][6 - 1] = c1 * ap3;
	   if (ndm < 3) return;

		   //!Compute jacobian transformation

		   //do j = 1, 3
			for(int j = 0;j< 3;j++){
				xs[j][1] = (xl[j][2] - xl[j][1])*shp[1][2]
					&          +(xl[j][3] - xl[j][4])*shp[1][3]
					&          +(xl[j][6] - xl[j][5])*shp[1][6]
					&          +(xl[j][7] - xl[j][8])*shp[1][7]

					//ay
					xs[j][2] = (xl[j][3] - xl[j][2])*shp[2][3]
					&          +(xl[j][4] - xl[j][1])*shp[2][4]
					&          +(xl[j][7] - xl[j][6])*shp[2][7]
					&          +(xl[j][8] - xl[j][5])*shp[2][8]

					xs[j][3] = (xl[j][5] - xl[j][1])*shp[3][5]
					&          +(xl[j][6] - xl[j][2])*shp[3][6]
					&          +(xl[j][7] - xl[j][3])*shp[3][7]
					&          +(xl[j][8] - xl[j][4])*shp[3][8]
					//end do
			}


			   //!Compute adjoint to jacobian

			   ad(1, 1) = xs(2, 2)*xs(3, 3) - xs(2, 3)*xs(3, 2)
			   ad(1, 2) = xs(3, 2)*xs(1, 3) - xs(3, 3)*xs(1, 2)
			   ad(1, 3) = xs(1, 2)*xs(2, 3) - xs(1, 3)*xs(2, 2)

			   ad(2, 1) = xs(2, 3)*xs(3, 1) - xs(2, 1)*xs(3, 3)
			   ad(2, 2) = xs(3, 3)*xs(1, 1) - xs(3, 1)*xs(1, 3)
			   ad(2, 3) = xs(1, 3)*xs(2, 1) - xs(1, 1)*xs(2, 3)

			   ad(3, 1) = xs(2, 1)*xs(3, 2) - xs(2, 2)*xs(3, 1)
			   ad(3, 2) = xs(3, 1)*xs(1, 2) - xs(3, 2)*xs(1, 1)
			   ad(3, 3) = xs(1, 1)*xs(2, 2) - xs(1, 2)*xs(2, 1)

			   //!Compute determinant of jacobian

			   xsj = xs(1, 1)*ad(1, 1) + xs(1, 2)*ad(2, 1) + xs(1, 3)*ad(3, 1)
			   rxsj = 1.d0 / xsj

			   //!Compute jacobian inverse

			   do j = 1, 3
				   do i = 1, 3
					   xs(i, j) = ad(i, j)*rxsj
					   end do
					   end do

					   //!Compute derivatives with repect to global coords.

					   do k = 1, 8

						   c1 = shp(1, k)*xs(1, 1) + shp(2, k)*xs(2, 1) + shp(3, k)*xs(3, 1)
						   c2 = shp(1, k)*xs(1, 2) + shp(2, k)*xs(2, 2) + shp(3, k)*xs(3, 2)
						   c3 = shp(1, k)*xs(1, 3) + shp(2, k)*xs(2, 3) + shp(3, k)*xs(3, 3)

						   shp(1, k) = c1
						   shp(2, k) = c2
						   shp(3, k) = c3

						   end do

						   end


*/
   return;
}
