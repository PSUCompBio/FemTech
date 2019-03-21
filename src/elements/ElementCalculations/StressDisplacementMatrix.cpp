#include "FemTech.h"
#include "blas.h"

void StressDisplacementMatrix(int e,int gp){

	//FT = F.transpose();
	if(ndim == 2){
		// B in 2D is a 3 x 2 matrix
		double B11,B12,B21,B22,B31,B32;
		double dnIdx = 0.0;
		double dnIdy = 0.0;
		double F11 = 0.0;
		double F21 = 0.0;
		double F12 = 0.0;
		double F22 = 0.0;

		B11 = dnIdx * F11;
		B12 = dnIdx * F21;

		B21 = dnIdy * F12;
		B22 = dnIdy * F22;

		B31 = dnIdx * F12 + dnIdy * F11;
		B32 = dnIdx * F22 + dnIdy * F21;
	}
	if(ndim == 3){
		// B matrix in 3D is a 6 x 3 matrix
		double B11,B12,B13,B21,B22,B23,B31,B32,B33;
		double dnIdx = 0.0;
		double dnIdy = 0.0;
		double dnIdz = 0.0;

		double F11 = 0.0;
		double F21 = 0.0;
		double F31 = 0.0;

		double F12 = 0.0;
		double F22 = 0.0;
		double F32 = 0.0;

		double F13 = 0.0;
		double F23 = 0.0;
		double F33 = 0.0;

		B11 = dnIdx * F11;
		B12 = dnIdx * F21;
		B13 = dnIdx * F31;

		B21 = dnIdy * F12;
		B22 = dnIdy * F22;
		B23 = dnIdy * F32;

		B31 = dnIdy*F13+dnIdz*F12 + dnIdy*F21+dnIdz*F22 + dnIdy*F33+dnIdz*F32;
		B32 = dnIdx*F13+dnIdz*F11 + dnIdx*F23+dnIdz*F21 + dnIdx*F33+dnIdz*F31;
		B33 = dnIdx*F12+dnIdy*F11 + dnIdx*F22+dnIdy*F21 + dnIdx*F32+dnIdy*F31;
	}

	return;
}
