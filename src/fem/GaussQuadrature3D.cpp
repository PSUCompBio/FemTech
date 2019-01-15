#include "digitalbrain.h"
void GaussQuadrature3D(int QuadratureRule, double *Chi,double *GaussWeights){
	/*-[--.----+----.----+----.-----------------------------------------]
	!      Purpose: Gauss quadrature for 3-d element

	!      Inputs:
	!         ll     - Number of points/direction

	!      Outputs:
	!         lint   - Total number of quadrature points
	!         s(4,*) - Gauss points (1-3) and weights (4)
	!-----[--.----+----.----+----.-----------------------------------------*/
   //printf("Rank %d: Guass quadrature for 3D elements!!\n", world_rank);


   //  2 x 2 x 2 pt. quadrature
	 if(QuadratureRule == 2){
		/*integration point 1*/
		 Chi[0*ndim + 0] = -0.577350269189626;
		 Chi[0*ndim + 1] = -0.577350269189626;
		 Chi[0*ndim + 2] = -0.577350269189626;
		/*integration point  2*/
		 Chi[1 * ndim + 0] = 0.577350269189626;
		 Chi[1 * ndim + 1] = -0.577350269189626;
		 Chi[1 * ndim + 2] = -0.577350269189626;
		/*integration point  3*/
		 Chi[2 * ndim + 0] = 0.577350269189626;
		 Chi[2 * ndim + 1] = 0.577350269189626;
		 Chi[2 * ndim + 2] = -0.577350269189626;
		/*integration point  4*/
		 Chi[3 * ndim + 0] = -0.577350269189626;
		 Chi[3 * ndim + 1] = 0.577350269189626;
		 Chi[3 * ndim + 2] = -0.577350269189626;
		/*integration point  5*/
		 Chi[4 * ndim + 0] = -0.577350269189626;
		 Chi[4 * ndim + 1] = -0.577350269189626;
		 Chi[4 * ndim + 2] = 0.577350269189626;
		/*integration point  6*/
		 Chi[5 * ndim + 0] = 0.577350269189626;
		 Chi[5 * ndim + 1] = -0.577350269189626;
		 Chi[5 * ndim + 2] = 0.577350269189626;
		/*integration point  7*/
		 Chi[6 * ndim + 0] = 0.577350269189626;
		 Chi[6 * ndim + 1] = 0.577350269189626;
		 Chi[6 * ndim + 2] = 0.577350269189626;
		/*integration point  8*/
		 Chi[7 * ndim + 0] = -0.577350269189626;
		 Chi[7 * ndim + 1] = 0.577350269189626;
		 Chi[7 * ndim + 2] = 0.577350269189626;

		 GaussWeights[0] = 1;
		 GaussWeights[1] = 1;
		 GaussWeights[2] = 1;
		 GaussWeights[3] = 1;
		 GaussWeights[4] = 1;
		 GaussWeights[5] = 1;
		 GaussWeights[6] = 1;
		 GaussWeights[7] = 1;
		 
	 }

   return;
}
