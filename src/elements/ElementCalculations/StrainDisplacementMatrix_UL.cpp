#include "FemTech.h"

void StrainDisplacementMatrix(int e, int gp, int nI, double *Bmat) {

	//FT = F.transpose();
	if(ndim == 2) {
		// Bmat in 2D is a 3 x 2 matrix
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
    const int indexStart = dsptr[e]+(gp*GaussPoints[e]+nI)*ndim;

		double dnIdx = dshp[indexStart];
		double dnIdy = dshp[indexStart+1];
		double dnIdz = dshp[indexStart+2];

		double B11,B12,B13,B21,B22,B23,B31,B32,B33;
		double B41,B42,B43,B51,B52,B53,B61,B62,B63;
    // Following Belytschko pg. 214 B_I^0 equation : 4.9.25
		B11 = dnIdx;
		B12 = 0;
		B13 = 0;

		B21 = 0;
		B22 = dnIdy;
		B23 = 0;

		B31 = 0;
		B32 = 0;
		B33 = dnIdz;

		B41 = 0;
    B42 = dnIdz;
    B43 = dnIdy;

		B51 = dnIdz;
    B52 = 0;
    B53 = dnIdx;

		B61 = dnIdy;
    B62 = dnIdx;
    B63 = 0;

    Bmat[0] = B11; Bmat[1] = B21; Bmat[2] = B31; Bmat[3] = B41; Bmat[4] = B51; Bmat[5] = B61;
    Bmat[6] = B12; Bmat[7] = B22; Bmat[8] = B32; Bmat[9] = B42; Bmat[10] = B52; Bmat[11] = B62;
    Bmat[12] = B13; Bmat[13] = B23; Bmat[14] = B33; Bmat[15] = B43; Bmat[16] = B53; Bmat[17] = B63;
	}
	return;
}
