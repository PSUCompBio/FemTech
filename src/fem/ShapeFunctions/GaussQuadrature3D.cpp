#include "FemTech.h"
void GaussQuadrature3D(int element, int nGaussPoint, double *Chi,double *GaussWeights){
	/*-[--.----+----.----+----.-----------------------------------------]
	!      Purpose: Gauss quadrature for 3-d element

	!      Inputs:
	!         ll     - Number of points/direction

	!      Outputs:
	!         lint   - Total number of quadrature points
	!         s(4,*) - Gauss points (1-3) and weights (4)
	!-----[--.----+----.----+----.-----------------------------------------*/
   //printf("Rank %d: Guass quadrature for 3D elements!!\n", world_rank);

	// 8-noded hex with 8 integration points
   //  2 x 2 x 2 pt. quadrature
	if (strcmp(ElementType[element], "C3D8") == 0 && nGaussPoint==8) {
#if 0
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
#endif
		 /*integration point 1*/
		 Chi[0 * ndim + 0] = -0.577350269189626;
		 Chi[0 * ndim + 1] = -0.577350269189626;
		 Chi[0 * ndim + 2] = 0.577350269189626;
		 /*integration point  2*/
		 Chi[1 * ndim + 0] = 0.577350269189626;
		 Chi[1 * ndim + 1] = -0.577350269189626;
		 Chi[1 * ndim + 2] = 0.577350269189626;
		 /*integration point  3*/
		 Chi[2 * ndim + 0] = 0.577350269189626;
		 Chi[2 * ndim + 1] = 0.577350269189626;
		 Chi[2 * ndim + 2] = 0.577350269189626;
		 /*integration point  4*/
		 Chi[3 * ndim + 0] = -0.577350269189626;
		 Chi[3 * ndim + 1] = 0.577350269189626;
		 Chi[3 * ndim + 2] = 0.577350269189626;
		 /*integration point  5*/
		 Chi[4 * ndim + 0] = -0.577350269189626;
		 Chi[4 * ndim + 1] = -0.577350269189626;
		 Chi[4 * ndim + 2] = -0.577350269189626;
		 /*integration point  6*/
		 Chi[5 * ndim + 0] = 0.577350269189626;
		 Chi[5 * ndim + 1] = -0.577350269189626;
		 Chi[5 * ndim + 2] = -0.577350269189626;
		 /*integration point  7*/
		 Chi[6 * ndim + 0] = 0.577350269189626;
		 Chi[6 * ndim + 1] = 0.577350269189626;
		 Chi[6 * ndim + 2] = -0.577350269189626;
		 /*integration point  8*/
		 Chi[7 * ndim + 0] = -0.577350269189626;
		 Chi[7 * ndim + 1] = 0.577350269189626;
		 Chi[7 * ndim + 2] = -0.577350269189626;

		 GaussWeights[0] = 1.0;
		 GaussWeights[1] = 1.0;
		 GaussWeights[2] = 1.0;
		 GaussWeights[3] = 1.0;
		 GaussWeights[4] = 1.0;
		 GaussWeights[5] = 1.0;
		 GaussWeights[6] = 1.0;
		 GaussWeights[7] = 1.0;
		 
	}

	// 4-noded tet with 1 integration point
	if (strcmp(ElementType[element], "C3D4") == 0 && nGaussPoint==1) {
		/*integration point 1*/
		Chi[0 * ndim + 0] = 0.25;
		Chi[0 * ndim + 1] = 0.25;
		Chi[0 * ndim + 2] = 0.25;

		GaussWeights[0] = 1.0/6.0;
	}

	// 4-noded tet with 4 integration point
	if (strcmp(ElementType[element], "C3D4") == 0 && nGaussPoint == 4) {
		/*integration point 1*/
		Chi[0 * ndim + 0] = 0.25;
		Chi[0 * ndim + 1] = 0.25;
		Chi[0 * ndim + 2] = 0.25;
		/*integration point 2*/
		Chi[1 * ndim + 0] = 0.25;
		Chi[1 * ndim + 1] = 0.25;
		Chi[1 * ndim + 2] = 0.25;
		/*integration point 3*/
		Chi[2 * ndim + 0] = 0.25;
		Chi[2 * ndim + 1] = 0.25;
		Chi[2 * ndim + 2] = 0.25;
		/*integration point 4*/
		Chi[3 * ndim + 0] = 0.25;
		Chi[3 * ndim + 1] = 0.25;
		Chi[3 * ndim + 2] = 0.25;

		GaussWeights[0] = 1.0 / 6.0;
		GaussWeights[1] = 1.0 / 6.0;
		GaussWeights[2] = 1.0 / 6.0;
		GaussWeights[3] = 1.0 / 6.0;
	}

   return;
}
