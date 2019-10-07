#include "FemTech.h"
#include "blas.h"

void TrussStressForceUpdate(int e, int gp, double *force){

	int x = 0, y = 1, z = 2;
	int node1 = 0, node2 = 1;
	double x01,x02,y01,y02,z01,z02;
	double x1,x2,y1,y2,z1,z2;
	int index=eptr[e];
	x01 = coordinates[ndim*connectivity[index+node1]+x];
	x02 = coordinates[ndim*connectivity[index+node2]+x];
	y01 = coordinates[ndim*connectivity[index+node1]+y];
	y02 = coordinates[ndim*connectivity[index+node2]+y];
	z01 = coordinates[ndim*connectivity[index+node1]+z];
	z02 = coordinates[ndim*connectivity[index+node2]+z];

	x1 = coordinates[ndim*connectivity[index+node1]+x]+displacements[ndim*connectivity[index+node1]+x];
	x2 = coordinates[ndim*connectivity[index+node2]+x]+displacements[ndim*connectivity[index+node2]+x];
	y1 = coordinates[ndim*connectivity[index+node1]+y]+displacements[ndim*connectivity[index+node1]+y];
	y2 = coordinates[ndim*connectivity[index+node2]+y]+displacements[ndim*connectivity[index+node2]+y];
	z1 = coordinates[ndim*connectivity[index+node1]+z]+displacements[ndim*connectivity[index+node1]+z];
	z2 = coordinates[ndim*connectivity[index+node2]+z]+displacements[ndim*connectivity[index+node2]+z];
  double le0= sqrt( pow((x02-x01),2.0) + pow((y02-y01),2.0) + pow((z02-z01),2.0) );
	double le = sqrt( pow((x2-x1),  2.0) +  pow((y2-y1), 2.0) + pow((z2-z1),  2.0) );
	// axial stretch
	double lambda = le/le0;
	FILE_LOG_SINGLE(DEBUGLOG, "coords = %3.3e\n", coordinates[ndim*connectivity[index+node1]+x]);
	FILE_LOG_SINGLE(DEBUGLOG, "disp = %3.3e\n", displacements[ndim*connectivity[index+node1]+x]);
	FILE_LOG_SINGLE(DEBUGLOG, "x2 = %3.3f\n",x2 );
  FILE_LOG_SINGLE(DEBUGLOG, "x1 = %3.3f\n",x1 );
	FILE_LOG_SINGLE(DEBUGLOG, "(x2-x1)^2 = %3.3f\n",pow((x2-x1),  2.0) );
	FILE_LOG_SINGLE(DEBUGLOG, "(y2-y1)^2 = %3.3f\n",pow((y2-y1), 2.0));
  FILE_LOG_SINGLE(DEBUGLOG, "(z2-z1)^2 = %3.3f\n",pow((z2-z1), 2.0));
  FILE_LOG_SINGLE(DEBUGLOG, "le0 = %3.3f\n",le0);
  FILE_LOG_SINGLE(DEBUGLOG, "le = %3.3f\n",le);
	FILE_LOG_SINGLE(DEBUGLOG, "lambda = %3.3f\n",lambda);
	// direction cosines
	double l = (x2-x1)/le;
	double m = (y2-y1)/le;
	double n = (z2-z1)/le;

	FILE_LOG_SINGLE(DEBUGLOG, "l = %3.3f\n",l);
	FILE_LOG_SINGLE(DEBUGLOG, "m = %3.3f\n",m);
  FILE_LOG_SINGLE(DEBUGLOG, "n = %3.3f\n",n);
	//create transformation matrix. following logan, section 3.7, page 93, equation 3.7.6
	double T[6];
	memset(T, 0.0, sizeof(T));
	T[0]=-1.0*l;
	T[1]=-1.0*m;
	T[2]=-1.0*n;
	T[3]=l;
	T[4]=m;
	T[5]=n;
	//T[0][1]=0;T[1][1]=0;T[2][1]=0;T[3][1]=l;T[4][1]=m;T[5][1]=n;

	//create local displacemtn vector
	double d[6];
	memset(d, 0.0, sizeof(d));
	d[0]=displacements[ndim*connectivity[index+node1]+x];
	d[1]=displacements[ndim*connectivity[index+node1]+y];
	d[2]=displacements[ndim*connectivity[index+node1]+z];
	d[3]=displacements[ndim*connectivity[index+node2]+x];
	d[4]=displacements[ndim*connectivity[index+node2]+y];
	d[5]=displacements[ndim*connectivity[index+node2]+z];

	// Muitiple T * d
	double sum = 0.0;
	for (int i = 0; i < 6; i++) {
		sum = sum + T[i] * d[i];
	}

  //green Lagrange strain
  double gls=0.5*(pow(lambda,2.0)-1.0);
	double area=1.0;
	double youngs = 1.0;
  double cauchystress = youngs * gls;
	double truss_axial_force  = cauchystress*area;
	FILE_LOG_SINGLE(DEBUGLOG, "gls = %3.3f\n",gls);

	FILE_LOG_SINGLE(DEBUGLOG, "s = %3.3f, f = %3.3f\n",cauchystress, truss_axial_force);

	// 6 values saved per gauss point for 3d
	// in voigt notation, sigma11
	pk2[pk2ptr[e] + 6 * gp + 0] = cauchystress;
	// in voigt notation, sigma22
	pk2[pk2ptr[e] + 6 * gp + 1] = 0;
	// in voigt notation, sigma33
	pk2[pk2ptr[e] + 6 * gp + 2] = 0;
	// in voigt notation, sigma23
	pk2[pk2ptr[e] + 6 * gp + 3] = 0;
	// in voigt notation, sigma13
	pk2[pk2ptr[e] + 6 * gp + 4] = 0;
	// in voigt notation, sigma12
	pk2[pk2ptr[e] + 6 * gp + 5] = 0;

	force[0] = truss_axial_force;
	force[1] = 0.0;
	force[2] = 0.0;
	force[3] = 0.0;
	force[4] = 0.0;
	force[5] = 0.0;
	return ;
}
