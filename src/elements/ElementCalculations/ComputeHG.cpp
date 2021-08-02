#include "FemTech.h"

#include "blas.h"

void ComputeHG(int e, double *force) {
	const int L[32] = {1, 1,-1,-1,-1,-1,1,1,1,-1,-1,1,-1,1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,-1,1,-1,1};
	double gamma[32];
	double temp[12];
	double HGforce[24];
	double ah;
	double Q = 0.1;
	int partid = pid[e];
	double rho = properties[MAXMATPARAMS * partid + 0];
	double ce = CalculateWaveSpeed(partid);
	double elementCoordinates[24];
	double v;
	double h[12];
	double g[12];
	for(int i = 0; i<12; i++){
		h[i]=0;
		g[i]=0;
		HGforce[i]=0;
		HGforce[i+12]=0;
		temp[i]=0;
		gamma[i]=0;
		gamma[i+12]=0;
		gamma[1+24]=0;
	}

	//Based on HG in LSDYNA
	for(int i = 0; i<32; i++){
		gamma[i]=0;
	}
	for (int i = eptr[e], j = 0; i < eptr[e+1]; ++i, ++j) {
		for (int k = 0; k < 3; ++k) {
			int index = ndim*connectivity[i]+k;
			elementCoordinates[j*ndim+k] = coordinates[index]+displacements[index];
		}
		//The following computation is for a different type of HG
	//	for(int alpha=0; alpha<4; alpha++){
	//	h[alpha*3+0] += velocities_half[ndim*connectivity[i]+0]*L[alpha*8+j];
	//	h[alpha*3+1] += velocities_half[ndim*connectivity[i]+1]*L[alpha*8+j];
	//	h[alpha*3+2] += velocities_half[ndim*connectivity[i]+2]*L[alpha*8+j];
	//}
	}
	{
	for(int i = eptr[e], j = 0; i < eptr[e+1]; ++i, ++j){
			int index = ndim*connectivity[i];
			for(int alpha=0; alpha<4; alpha++){
				temp[alpha*3+0] += (coordinates[index+0]+displacements[index+0])*L[alpha*8+j];
				temp[alpha*3+1] += (coordinates[index+1]+displacements[index+1])*L[alpha*8+j];
				temp[alpha*3+2] += (coordinates[index+2]+displacements[index+2])*L[alpha*8+j];
			}
		}
	const unsigned int indexStart = dsptr[e] + GaussPoints[e] * ndim;
	for(int i = eptr[e], j = 0; i < eptr[e+1]; ++i, ++j){
		const unsigned int indexShp = dsptr[e] + j * ndim;
		for(int alpha=0; alpha<4; alpha++){
			gamma[alpha*8+j] = L[alpha*8+j] - (dshp[indexShp]*F_XiInverse[0] + dshp[indexShp+1]*F_XiInverse[3] + dshp[indexShp+2]*F_XiInverse[6])*temp[alpha*3+0] - (dshp[indexShp]*F_XiInverse[1] +
		      dshp[indexShp+1]*F_XiInverse[4] + dshp[indexShp +2]*F_XiInverse[7])*temp[alpha*3+1] - (dshp[indexShp]*F_XiInverse[2] + dshp[indexShp+1]*F_XiInverse[5] + dshp[indexShp +2]*F_XiInverse[8])*temp[alpha*3+2];
			}
		}
	}
	for (int i = eptr[e], j = 0; i < eptr[e+1]; ++i, ++j)
		for(int alpha=0; alpha<4; alpha++){
		g[alpha*3+0] += velocities_half[ndim*connectivity[i]+0]*gamma[alpha*8+j];
		g[alpha*3+1] += velocities_half[ndim*connectivity[i]+1]*gamma[alpha*8+j];
		g[alpha*3+2] += velocities_half[ndim*connectivity[i]+2]*gamma[alpha*8+j];
		}
	v = volumeHexahedron(elementCoordinates);	
	ah = (Q*rho*pow(v,2.0/3.0)*ce)/4;
	for (int i = eptr[e], j = 0; i < eptr[e+1]; ++i, ++j) {
		for(int alpha=0; alpha<4; alpha++){
			HGforce[ndim*j+0] += ah*g[alpha*3+0]*gamma[alpha*8+j];
			HGforce[ndim*j+1] += ah*g[alpha*3+1]*gamma[alpha*8+j];
			HGforce[ndim*j+2] += ah*g[alpha*3+2]*gamma[alpha*8+j];
			}
	}
	for (int i = eptr[e], j = 0; i < eptr[e+1]; ++i, ++j) {
		force[ndim*j+0] += HGforce[ndim*j+0];
		force[ndim*j+1] += HGforce[ndim*j+1];
		force[ndim*j+2] += HGforce[ndim*j+2];
		f_hg[connectivity[i+j]*ndim + 0] += HGforce[ndim*j+0];
		f_hg[connectivity[i+j]*ndim + 1] += HGforce[ndim*j+1];
		f_hg[connectivity[i+j]*ndim + 2] += HGforce[ndim*j+2];	
		}
}
