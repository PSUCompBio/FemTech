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
	sigma_n[sigmaptr[e] + 6 * gp + 0] = cauchystress;
	// in voigt notation, sigma22
	sigma_n[sigmaptr[e] + 6 * gp + 1] = 0;
	// in voigt notation, sigma33
	sigma_n[sigmaptr[e] + 6 * gp + 2] = 0;
	// in voigt notation, sigma23
	sigma_n[sigmaptr[e] + 6 * gp + 3] = 0;
	// in voigt notation, sigma13
	sigma_n[sigmaptr[e] + 6 * gp + 4] = 0;
	// in voigt notation, sigma12
	sigma_n[sigmaptr[e] + 6 * gp + 5] = 0;

	force[0] = truss_axial_force;
	force[1] = 0.0;
	force[2] = 0.0;
	force[3] = 0.0;
	force[4] = 0.0;
	force[5] = 0.0;
	//
  // // Calculate Tranformation
	// double dir_truss_x[3];
	// double dir_truss_y[3];
	// double dir_truss_z[3];
	// dir_truss_x[0] = x2 - x1;
	// dir_truss_x[1] = y2 - y1;
	// dir_truss_x[2] = z2 - z1;
	// if (dir_truss_x[0] != 0) {
	// 	dir_truss_y[1] = 1;
	// 	dir_truss_y[2] = 0;
	// 	dir_truss_y[0] = -(((dir_truss_y[1] * dir_truss_x[1]) + (dir_truss_y[2] * dir_truss_x[2])) / dir_truss_x[0]);
	// }
	// else {
	// 	if ((dir_truss_x[0] == 0) && (dir_truss_x[1] != 0)) {
	// 		dir_truss_y[0] = 1;
	// 		dir_truss_y[2] = 1;
	// 		dir_truss_y[1] = -(((dir_truss_y[2] * dir_truss_x[2]) + (dir_truss_y[0] * dir_truss_x[0])) / dir_truss_x[1]);
	// 	}
	// 	else if ((dir_truss_x[0] == 0) && (dir_truss_x[2] != 0)) {
	// 		dir_truss_y[0] = 0;
	// 		dir_truss_y[1] = 1;
	// 		dir_truss_y[2] = -(((dir_truss_y[0] * dir_truss_x[0]) + (dir_truss_y[1] * dir_truss_x[1])) / dir_truss_x[2]);
	// 	}
	// 	else {
	// 		printf("Truss element is not practically possible.\n");
	// 	}
	// }
  // // cross product: of x cross y to get z
	// double a1 = dir_truss_x[0];
	// double a2 = dir_truss_x[1];
	// double a3 = dir_truss_x[2];
	// double b1 = dir_truss_y[0];
	// double b2 = dir_truss_y[1];
	// double b3 = dir_truss_y[2];
  // dir_truss_z[0] = (a2*b3 - a3*b2);
	// dir_truss_z[1] = (a1*b3 - a3*b1);
	// dir_truss_z[2] = (a1*b1 - a2*b1);
	//
	// //normalize each vector
	// double norm_x = sqrt(dir_truss_x[0]*dir_truss_x[0] + dir_truss_x[1]*dir_truss_x[1] + dir_truss_x[2]*dir_truss_x[2] );
	// double norm_y = sqrt(dir_truss_y[0]*dir_truss_y[0] + dir_truss_y[1]*dir_truss_y[1] + dir_truss_y[2]*dir_truss_y[2] );
	// double norm_z = sqrt(dir_truss_z[0]*dir_truss_z[0] + dir_truss_z[1]*dir_truss_z[1] + dir_truss_z[2]*dir_truss_z[2] );
	// // x
  // dir_truss_x[0] =  dir_truss_x[0]/norm_x;
	// dir_truss_x[1] =  dir_truss_x[1]/norm_x;
	// dir_truss_x[2] =  dir_truss_x[2]/norm_x;
 	// // y
	// dir_truss_y[0] =  dir_truss_y[0]/norm_y;
	// dir_truss_y[1] =  dir_truss_y[1]/norm_y;
	// dir_truss_y[2] =  dir_truss_y[2]/norm_y;
	// // z
	// dir_truss_z[0] =  dir_truss_z[0]/norm_z;
	// dir_truss_z[1] =  dir_truss_z[1]/norm_z;
	// dir_truss_z[2] =  dir_truss_z[2]/norm_z;
	//
	// // define global coordinate system
	// double dir_global_x[3];
	// dir_global_x[0]=1;
	// dir_global_x[1]=0;
	// dir_global_x[2]=0;
	// double dir_global_y[3];
	// dir_global_y[0]=0;
	// dir_global_y[1]=1;
	// dir_global_y[2]=0;
	// double dir_global_z[3];
	// dir_global_z[0]=0;
	// dir_global_z[1]=0;
	// dir_global_z[2]=1;
	// //normalize each global vector - not really needed but left in b/c it can
	// // accomadate different global orientations in future  - IS THIS NEEDED??
	// double gnorm_x = sqrt(dir_global_x[0]*dir_global_x[0] + dir_global_x[1]*dir_global_x[1] + dir_global_x[2]*dir_global_x[2] );
	// double gnorm_y = sqrt(dir_global_y[0]*dir_global_y[0] + dir_global_y[1]*dir_global_y[1] + dir_global_y[2]*dir_global_y[2] );
	// double gnorm_z = sqrt(dir_global_z[0]*dir_global_z[0] + dir_global_z[1]*dir_global_z[1] + dir_global_z[2]*dir_global_z[2] );
	//
	// // define direction cosines
	// double dot_global_x_truss_x = dir_global_x[0]*dir_truss_x[0] + dir_global_x[1]*dir_truss_x[1] + dir_global_x[2]*dir_truss_x[2];
	// double dot_global_x_truss_y = dir_global_x[0]*dir_truss_y[0] + dir_global_x[1]*dir_truss_y[1] + dir_global_x[2]*dir_truss_y[2];
	// double dot_global_x_truss_z = dir_global_x[0]*dir_truss_z[0] + dir_global_x[1]*dir_truss_z[1] + dir_global_x[2]*dir_truss_z[2];
	// double dot_global_y_truss_x = dir_global_y[0]*dir_truss_x[0] + dir_global_y[1]*dir_truss_x[1] + dir_global_y[2]*dir_truss_x[2];
	// double dot_global_y_truss_y = dir_global_y[0]*dir_truss_y[0] + dir_global_y[1]*dir_truss_y[1] + dir_global_y[2]*dir_truss_y[2];
	// double dot_global_y_truss_z = dir_global_y[0]*dir_truss_z[0] + dir_global_y[1]*dir_truss_z[1] + dir_global_y[2]*dir_truss_z[2];
	// double dot_global_z_truss_x = dir_global_z[0]*dir_truss_x[0] + dir_global_z[1]*dir_truss_x[1] + dir_global_z[2]*dir_truss_x[2];
	// double dot_global_z_truss_y = dir_global_z[0]*dir_truss_y[0] + dir_global_z[1]*dir_truss_y[1] + dir_global_z[2]*dir_truss_y[2];
	// double dot_global_z_truss_z = dir_global_z[0]*dir_truss_z[0] + dir_global_z[1]*dir_truss_z[1] + dir_global_z[2]*dir_truss_z[2];
	// // l
	// double l_1 = dot_global_x_truss_x / (gnorm_x * norm_x);
	// double l_2 = dot_global_x_truss_y / (gnorm_x * norm_y);
	// double l_3 = dot_global_x_truss_z / (gnorm_x * norm_z);
  // // m
	// double m_1 = dot_global_y_truss_x / (gnorm_y * norm_x);
	// double m_2 = dot_global_y_truss_y / (gnorm_y * norm_y);
	// double m_3 = dot_global_y_truss_z / (gnorm_y * norm_z);
	// // n
	// double n_1 = dot_global_z_truss_x / (gnorm_z * norm_x);
	// double n_2 = dot_global_z_truss_y / (gnorm_z * norm_y);
	// double n_3 = dot_global_z_truss_z / (gnorm_z * norm_z);
	//
	// // choice - 1: Stress Transformation Matrix
	// // choice - 2: Strain Transformation Matrix
	// // choice - 3: 3X3 transformation matrix as produced in my thesis in electrophysiology chapter
	// // int choice = 3;
	// // if (choice == 2) {
	// // 	double T[6][6];
	// // 	memset(T, 0.0, sizeof(T));
	// // 	T[0][0]=l_1*l_1;T[1][0]=m_1*m_1;T[2][0]=n_1*n_1;T[3][0]=l_1*m_1;T[4][0]=m_1*n_1;T[5][0]=n_1*l_1;
	// // 	T[0][1]=l_2*l_2;T[1][1]=m_2*m_2;T[2][1]=n_2*n_2;T[3][1]=l_2*m_2;T[4][1]=m_2*n_2;T[5][1]=n_2*l_2;
	// // 	T[0][2]=l_3*l_3;T[1][2]=m_3*m_3;T[2][2]=n_3*n_3;T[3][2]=l_3*m_3;T[4][2]=m_3*n_3;T[5][2]=n_3*l_3;
	// // 	T[0][3]=2.0*l_1*l_2;T[1][3]=2.0*m_1*m_2;T[2][3]=2.0*n_1*n_2;T[3][3]=(l_1*m_2)+(l_2*m_1);T[4][3]=(m_1*n_2)+(m_2*n_1);T[5][3]=(l_1*n_2)+(l_2*n_1);
	// // 	T[0][4]=2.0*l_2*l_3;T[1][4]=2.0*m_2*m_3;T[2][4]=2.0*n_2*n_3;T[3][4]=(l_2*m_3)+(l_3*m_2);T[4][4]=(m_2*n_3)+(m_3*n_2);T[5][4]=(l_2*n_3)+(l_3*n_2);
	// // 	T[0][5]=2.0*l_1*l_3;T[1][5]=2.0*m_1*m_3;T[2][5]=2.0*n_1*n_3;T[3][5]=(l_1*m_3)+(l_3*m_1);T[4][5]=(m_1*n_3)+(m_3*n_1);T[5][5]=(l_1*n_3)+(l_3*n_1);
	// 	// T.row(0) << (l_1 * l_1), (m_1 * m_1), (n_1 * n_1), (l_1 * m_1), (m_1 * n_1), (n_1 * l_1);
	// 	// T.row(1) << (l_2 * l_2), (m_2 * m_2), (n_2 * n_2), (l_2 * m_2), (m_2 * n_2), (n_2 * l_2);
	// 	// T.row(2) << (l_3 * l_3), (m_3 * m_3), (n_3 * n_3), (l_3 * m_3), (m_3 * n_3), (n_3 * l_3);
	// 	// T.row(3) << (2 * l_1 * l_2), (2 * m_1 * m_2), (2 * n_1 * n_2), ((l_1 * m_2) + (l_2 * m_1)), ((m_1 * n_2) + (m_2 * n_1)), ((l_1 * n_2) + (l_2 * n_1));
	// 	// T.row(4) << (2 * l_2 * l_3), (2 * m_2 * m_3), (2 * n_2 * n_3), ((l_2 * m_3) + (l_3 * m_2)), ((m_2 * n_3) + (m_3 * n_2)), ((l_2 * n_3) + (l_3 * n_2));
	// 	// T.row(5) << (2 * l_1 * l_3), (2 * m_1 * m_3), (2 * n_1 * n_3), ((l_1 * m_3) + (l_3 * m_1)), ((m_1 * n_3) + (m_3 * n_1)), ((l_1 * n_3) + (l_3 * n_1));
	// // }
	// // else if (choice == 1) {
	// // 	double T[6][6];
	// // 	memset(T, 0.0, sizeof(T));
	// // 	T[0][0]=l_1*l_1;T[1][0]=m_1*m_1;T[2][0]=n_1*n_1; T[3][0]=2.0*l_1*m_1;T[4][0]=2.0*m_1*n_1;T[5][0]=2.0*n_1*l_1;
	// // 	T[0][1]=l_2*l_2;T[1][1]=m_2*m_2;T[2][1]=n_2*n_2; T[3][1]=2.0*l_2*m_2;T[4][1]=2.0*m_2*n_2;T[5][1]=2.0*n_2*l_2;
	// // 	T[0][2]=l_3*l_3;T[1][2]=m_3*m_3;T[2][2]=n_3*n_3; T[3][2]=2.0*l_3*m_3;T[4][2]=2.0*m_3*n_3;T[5][2]=2.0*n_3*l_3;
	// // 	T[0][3]=l_1*l_2;T[1][3]=m_1*m_2;T[2][3]=n_1*n_2;T[3][3]=(l_1*m_2)+(l_2*m_1);T[4][3]=(m_1*n_2)+(m_2*n_1);T[5][3]=(l_1*n_2)+(l_2*n_1);
	// // 	T[0][4]=l_2*l_3;T[1][4]=m_2*m_3;T[2][4]=n_2*n_3;T[3][4]=(l_2*m_3)+(l_3*m_2);T[4][4]=(m_2*n_3)+(m_3*n_2);T[5][4]=(l_2*n_3)+(l_3*n_2);
	// // 	T[0][5]=l_1*l_3;T[1][5]=m_1*m_3;T[2][5]=n_1*n_3;T[3][5]=(l_1*m_3)+(l_3*m_1);T[4][5]=(m_1*n_3)+(m_3*n_1);T[5][5]=(l_1*n_3)+(l_3*n_1);
	// // 	// T.row(0) << pow(l_1, 2), pow(m_1, 2), pow(n_1, 2), (2 * l_1 * m_1), (2 * m_1 * n_1), (2 * n_1 * l_1);
	// // 	// T.row(1) << pow(l_2, 2), pow(m_2, 2), pow(n_2, 2), (2 * l_2 * m_2), (2 * m_2 * n_2), (2 * n_2 * l_2);
	// // 	// T.row(2) << pow(l_3, 2), pow(m_3, 2), pow(n_3, 2), (2 * l_3 * m_3), (2 * m_3 * n_3), (2 * n_3 * l_3);
	// // 	// T.row(3) << (1 * l_1 * l_2), (1 * m_1 * m_2), (1 * n_1 * n_2), ((l_1 * m_2) + (l_2 * m_1)), ((m_1 * n_2) + (m_2 * n_1)), ((l_1 * n_2) + (l_2 * n_1));
	// // 	// T.row(4) << (1 * l_2 * l_3), (1 * m_2 * m_3), (1 * n_2 * n_3), ((l_2 * m_3) + (l_3 * m_2)), ((m_2 * n_3) + (m_3 * n_2)), ((l_2 * n_3) + (l_3 * n_2));
	// // 	// T.row(5) << (1 * l_1 * l_3), (1 * m_1 * m_3), (1 * n_1 * n_3), ((l_1 * m_3) + (l_3 * m_1)), ((m_1 * n_3) + (m_3 * n_1)), ((l_1 * n_3) + (l_3 * n_1));
	// // }
	// // else if (choice == 3) {
	// 	double T[3][3];
	// 	memset(T, 0.0, sizeof(T));
	// 	T[0][0]=l_1;T[1][0]=l_2;T[2][0]=l_3;
	// 	T[0][1]=m_1;T[1][1]=m_2;T[2][1]=m_3;
	// 	T[0][2]=n_1;T[1][2]=n_2;T[2][2]=n_3;
	// 	// T.row(0) << l_1, l_2, l_3;
	// 	// T.row(1) << m_1, m_2, m_3;
	// 	// T.row(2) << n_1, n_2, n_3;
	// //}
	// // else {
	// // 	printf("WRONG CHOICE OF TRANSFORMATION MATRIX\n");
	// // }
	//
	// // Compute inverse of T
	// double T_inv[3][3];
	// double T_inv_transpose[3][3];
	// memset(T, 0.0, sizeof(T_inv));
	// // first compute determinant
	// double tmp1 = T[0][0]*(T[1][1]*T[2][2] - T[1][2]*T[2][1]);
	// double tmp2 = T[0][1]*(T[1][0]*T[2][2] - T[1][2]*T[2][0]);
	// double tmp3 = T[0][2]*(T[1][0]*T[2][1] - T[1][1]*T[2][0]);
	// double detT = tmp1 - tmp2 + tmp3;
	// // now compute inverse
	// T_inv[0][0] = 1/detT*(T[1][1]*T[2][2] - T[1][2]*T[2][1]);
	// T_inv[0][1] = 1/detT*(T[0][2]*T[2][1] - T[0][1]*T[2][2]);
	// T_inv[0][2] = 1/detT*(T[0][1]*T[1][2] - T[0][2]*T[1][1]);
	//
	// T_inv[1][0] = 1/detT*(T[1][2]*T[2][0] - T[1][0]*T[2][2]);
	// T_inv[1][1] = 1/detT*(T[0][0]*T[2][2] - T[0][2]*T[2][0]);
	// T_inv[1][2] = 1/detT*(T[0][2]*T[1][0] - T[0][0]*T[1][2]);
	//
	// T_inv[2][0] = 1/detT*(T[1][0]*T[2][1] - T[1][1]*T[2][0]);
	// T_inv[2][1] = 1/detT*(T[0][1]*T[2][0] - T[0][0]*T[2][1]);
	// T_inv[2][2] = 1/detT*(T[0][0]*T[1][1] - T[0][1]*T[1][0]);
	// // Finding the transpose of matrix T_inv
  // for(int i=0; i<ndim; ++i){
  //   for(int j=0; j<ndim; ++j){
  //    	T_inv_transpose[j][i] = T_inv[i][j];
  //   }
	// }
	//
	// // Compute deformation gradient of truss  - assumes incompressibilty
	// // because F is typically stored at gauss points, the loop over gauss points
	// // is included. The size of F is defined in Shapefunctions.cpp. It is already
  // // created but it needs to referenced correctly.
	// index = fptr[e] + ndim*ndim*gp;
	// F[index+ndim*0+0] = lambda;
	// F[index+ndim*1+1] = 1.0/sqrt(lambda);
	// F[index+ndim*2+2] = 1.0/sqrt(lambda);
	//
	// // Muitiple T_inv_transpose * F
	// double sum = 0;
	// double multiply[3][3];
	// for (int c = 0; c < ndim; c++) {
	// 	for (int d = 0; d < ndim; d++) {
	// 		for (int k = 0; k < ndim; k++) {
	// 			sum = sum + T_inv_transpose[c][k]*F[index+ndim*d+k];
	// 		}
	// 		multiply[c][d] = sum;
	// 		sum = 0;
	// 	}
	// }

	// // mulitipl F * T_inv
	// sum = 0;
	// for (int c = 0; c < ndim; c++) {
	// 	for (int d = 0; d < ndim; d++) {
	// 		for (int k = 0; k < ndim; k++) {
	// 			sum = sum + multiply[c][k]*T_inv[d][k];
	// 		}
	// 		F[c][d] = sum;
	// 		sum = 0;
	// 	}
	// }
	return ;
}
