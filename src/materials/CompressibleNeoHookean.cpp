#include "FemTech.h"
#include "blas.h"

// plane strain or three-dimensional compressible neo-Hookean
// Evaluates the Cauchy stress tensor

void CompressibleNeoHookean(int e, int gp){

//From Bonet and Wood - Flagshyp
//mu              = properties(2);
//lambda          = properties(3);
//J               = kinematics.J;
//b               = kinematics.b;
//Cauchy          = (mu/J)*(b - cons.I) + (lambda/J)*log(J)*cons.I;

	double mu = 1.0;
 	double lambda = 1.0;
	//double J = 1;
	//printf("element %d, gauss point %d\n",e,gp);

	if(ndim == 2){
		// 6 values saved per gauss point for 3d
		for(int i=0;i<3;i++){
			int index = cptr[e]+3*gp+i;
			//cauchy[index] = 0.0
		  //printf("cauchy[%d]\n",index);
		}
	}
	if(ndim == 3){
		// 6 values saved per gauss point for 3d
		for(int i=0;i<6;i++){
			int index = cptr[e]+6*gp+i;
			//cauchy[index] = 0.0
			//printf("cauchy[%d]\n",index);
		}
	}
	//printf("\n");
	return;


}
