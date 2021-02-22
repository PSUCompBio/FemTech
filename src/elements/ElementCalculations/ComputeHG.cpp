#include "FemTech.h"

#include "blas.h"

void ComputeHG(int e, double *force) {
	const int L[32] = {1, 1,-1,-1,-1,-1,1,1,1,-1,-1,1,-1,1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,-1,1,-1,1};
	double HGforce[24];
	double ah;
	double Q = 0.1;
	double rho = properties[MAXMATPARAMS * pid[e] + 0];
	double lambda = properties[MAXMATPARAMS * pid[e] + 2];
	double mu = properties[MAXMATPARAMS * pid[e] + 1];
	double nu = 0.5*lambda/(lambda+mu);
	double elementCoordinates[24];
	double v;
	double h[12];
	for(int i = 0; i<12; i++){
		h[i]=0;
		HGforce[i]=0;
		HGforce[i+12]=0;
	}
	for (int i = eptr[e], j = 0; i < eptr[e+1]; ++i, ++j) {
		for (int k = 0; k < 3; ++k) {
			int index = ndim*connectivity[i]+k;
			elementCoordinates[j*ndim+k] = coordinates[index]+displacements[index];
		}
		for(int alpha=0; alpha<4; alpha++){
		h[alpha*3+0] += velocities_half[ndim*connectivity[i]+0]*L[alpha*8+i];
		h[alpha*3+1] += velocities_half[ndim*connectivity[i]+1]*L[alpha*8+i];
		h[alpha*3+2] += velocities_half[ndim*connectivity[i]+2]*L[alpha*8+i];
//		printf("%f %f %f\n", h[alpha*3+0], h[alpha*3+1], h[alpha*3+2]);
	}
//	printf("%.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n", h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7], h[8], h[9], h[10], h[11]);
//	printf("%f %f %f\n", velocities_half[ndim*connectivity[i]+0], velocities_half[ndim*connectivity[i]+1], velocities_half[ndim*connectivity[i]+2]);
	}
//	printf("%.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n", h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7], h[8], h[9], h[10], h[11]);
	v = volumeHexahedron(elementCoordinates);
	double ce = sqrt(lambda*(1.0/nu-1.0)/rho);
	ah = (Q*rho*pow(v,2.0/3.0)*ce)/4;
	for (int i = eptr[e]; i < eptr[e+1]; ++i) {
			for(int alpha=0; alpha<4; alpha++){
						HGforce[ndim*connectivity[i]+0] += ah*h[alpha*3+0]*L[alpha*8+i];
						HGforce[ndim*connectivity[i]+1] += ah*h[alpha*3+1]*L[alpha*8+i];
						HGforce[ndim*connectivity[i]+2] += ah*h[alpha*3+2]*L[alpha*8+i];
					//	printf("%.10f %.10f %.10f\n", h[alpha*3+0], h[alpha*3+1], h[alpha*3+2]);
	}
	}
	for (int i = eptr[e]; i < eptr[e+1]; ++i) {
		printf("%.10f %.10f %.10f %f first\n", force[ndim*connectivity[i]+0], force[ndim*connectivity[i]+1], force[ndim*connectivity[i]+2], Time);
		force[ndim*connectivity[i]+0] += HGforce[ndim*connectivity[i]+0];
		force[ndim*connectivity[i]+1] += HGforce[ndim*connectivity[i]+1];
		force[ndim*connectivity[i]+2] += HGforce[ndim*connectivity[i]+2];
		printf("%.10f %.10f %.10f %f second\n", HGforce[ndim*connectivity[i]+0], HGforce[ndim*connectivity[i]+1], HGforce[ndim*connectivity[i]+2], Time);
	}
}
