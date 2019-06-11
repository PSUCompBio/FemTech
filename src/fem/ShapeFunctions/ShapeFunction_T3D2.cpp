#include "FemTech.h"

void get_cos(int e, double *cx, double *cy, double *cz){
      // int e :: elements' id
      // double cx,cy,cz :: truss direction
      // int eptr[e+1] :: number of nodes per element
      // For each element e, you obtain the involved nodes looping untin e+1.
      double l0 = {0};
      double x[2] = {0};
      double y[2] = {0};
      double z[2] = {0};
      
      for (int k=0; k < eptr[e+1]; k++){
            x[k] = coordinates[ndim*connectivity[eptr[e]+k]+0];
            y[k] = coordinates[ndim*connectivity[eptr[e]+k]+1];
            z[k] = coordinates[ndim*connectivity[eptr[e]+k]+2];
      }
//      Checking the coordinates of each node, for element id:e

       for (int k=0; k < eptr[e+1]; k++){
            printf("nodes per elem: %d, ide node: %d\n",eptr[e+1],k);
            printf("coordinates node: %d, are [%5.3e,%5.3e,%5.3e]\n",k,x[k],y[k],y[k]);
       }
	 printf("\n");
      
      l0 = sqrt(pow(x[1]-x[0],2) + pow(y[1]-y[0],2) + pow(z[1]-z[0],2));
      
      printf("length of the truss: %5.3e\n",l0);
      
      *cx = ((x[1] - x[0]) / l0);
      *cy = ((y[1] - y[0]) / l0);
      *cz = ((z[1] - z[0]) / l0);
      
      printf("the truss direction: [%5.3e,%5.3e,%5.3e]\n",*cx,*cy,*cz);
      
      return;
}

//https://www.sanfoundry.com/cpp-programming-examples-numerical-problems-algorithms/
//https://www.sanfoundry.com/c-programming-examples-numerical-problems-algorithms/
//https://www.sanfoundry.com/1000-c-algorithms-problems-programming-examples/

void ShapeFunction_T3D2(int e, int gp, double *Chi, double *detJ){
      // int e :: element ID
      // int gp :: Gauss point ID
      // double *Chi :: Gauss points array coordinates
      // double *detJ :: Jacs determinant array

	double chi, eta, iota;
	chi =  Chi[ndim*gp + 0]; // local x
	eta =  Chi[ndim*gp + 1]; // local y
	iota = Chi[ndim*gp + 2]; // local z

	// The shape functions
	int g = GaussPoints[e];
	const int indexShp = gptr[e] + gp * g;
      const int indexDshp = dsptr[e] + gp * g * ndim;

	shp[indexShp + 0] = 0.5*(1.0-chi);
	shp[indexShp + 1] = 0.5*(1.0+chi);

	// The first derivatives

	// with respect to local x
	dshp[indexDshp + ndim * 0 + 0] = -0.5; // derivative of N_1
 	dshp[indexDshp + ndim * 1 + 0] =  0.5; // derivative of N_2

	// with respect to local y
	dshp[indexDshp + ndim * 0 + 1] =  0.0; // derivative of N_1
	dshp[indexDshp + ndim * 1 + 1] =  0.0; // derivative of N_2
	
	// with respect to local z
	dshp[indexDshp + ndim * 0 + 2] =  0.0; // derivative of N_1
	dshp[indexDshp + ndim * 1 + 2] =  0.0; // derivative of N_2
	 
	// Here we start the transformation 6-dof to 2-dof.
	 
	double cx,cy,cz;
	double T[12];
	
	cx = 0.0;
	cy = 0.0;
	cz = 0.0;
	
	printf("the initial truss direction: [%5.3e,%5.3e,%5.3e]\n",cx,cy,cz);
	
	get_cos(e, &cx, &cy, &cz);
	
	printf("the truss direction: [%5.3e,%5.3e,%5.3e]\n",cx,cy,cz);
	
	printf("Number of spatial dimensions: %d\n",ndim);
	
	// Compute the Jacobian's determinant
	
	double xl[2] = {0};
      double x[2] = {0};
      double y[2] = {0};
      double z[2] = {0};
      
      for (int k=0; k < eptr[e+1]; k++){
            x[k] = coordinates[ndim*connectivity[eptr[e]+k]+0];
            y[k] = coordinates[ndim*connectivity[eptr[e]+k]+1];
            z[k] = coordinates[ndim*connectivity[eptr[e]+k]+2];
      }
	
	for (int k=0; k< eptr[e+1]; k++){
	      xl[k] = x[k]*cx + y[k]*cy + z[k]*cz;
	}
	
	printf("Elem id: %d, coord-0: %5.3e, coord-1: %5.3e\n",e,xl[0],xl[1]);
	
	detJ[gp] = xl[1] - xl[0];
	
	// The Global first derivatives. Only for direction x. 

	// with respect to local x
	dshp[indexDshp + ndim * 0 + 0] = dshp[indexDshp + ndim * 0 + 0] / detJ[gp]; // derivative of N_1
 	dshp[indexDshp + ndim * 1 + 0] = dshp[indexDshp + ndim * 1 + 0] / detJ[gp]; // derivative of N_2

	// //for debugging can be removed...
	// if (debug && 1==0) {
  //   printf("DEBUG J Inv\n");
	// 	printf("%8.4e %8.4e %8.4e\n", J_Inv[0], J_Inv[3], J_Inv[6]);
	// 	printf("%8.4e %8.4e %8.4e\n", J_Inv[1], J_Inv[4], J_Inv[7]);
	// 	printf("%8.4e %8.4e %8.4e\n", J_Inv[2], J_Inv[5], J_Inv[8]);
	// 	printf("\n");
	// }
   return;
}
