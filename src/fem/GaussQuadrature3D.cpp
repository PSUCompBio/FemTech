#include "digitalbrain.h"
double GaussQuadrature3D(int QuadratureRule, double []){
	/*-[--.----+----.----+----.-----------------------------------------]
	!      Purpose: Gauss quadrature for 3-d element

	!      Inputs:
	!         ll     - Number of points/direction

	!      Outputs:
	!         lint   - Total number of quadrature points
	!         s(4,*) - Gauss points (1-3) and weights (4)
	!-----[--.----+----.----+----.-----------------------------------------*/
   printf("Rank %d: Guass quadrature for 3D elements!!\n", world_rank);

	 double ig[4] = {-1, 1, 1, -1};
	 double jg[4] = {-1, -1, 1, 1};
	 //double sf[4];

   //  2 x 2 x 2 pt. quadrature
	 if(QuadratureRule == 2){
		 int lint = 8;
		 double g    = 1.0/sqrt(3.0);
		 for (int i=0;i<4;i++){
			 s[i] = 1.;
			 //sf[0][i]		 = ig[i]*g;
			 //sf[0][i+3]	 = sf[0][i];
			 //sf[1][i]		 = jg[i]*g;
			 //sf[1][i+3]	 = sf[1][i];
			 //sf[2][i]		 = g;
			 //sf[2][i+3]	 = -g;
			 //sf[3][i] 	 	= 1.0;
			 //sf[3][i+3]  = 1.0;
			 //printf("%2.5f %2.5f %2.5f %2.5f\n",sf[0][i],sf[1][i],sf[2][i],sf[3][i]);
		 }
	 }

   return;
}
