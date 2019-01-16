#include "FemTech.h"

/* Global Variables */
int *gptr;
int *GaussPoints;
double *shp;
int *nShapeFunctions;

void ShapeFunctions(){

	// Global Array - keeps track of how many gauss points there are 
	// per element.
	GaussPoints = (int *)malloc(nelements * sizeof(int));
	// Global Array - keeps track of how many shp functions there are 
	// per element.
	nShapeFunctions = (int *)malloc(nelements * sizeof(int));

	// Global array -
	// When number of shape functions in an element equals 
	// the number of nodes in an element, gptr array looks 
	// similar to eptr array. However, the number of quadrature
	// points (a.k.a gauss points) can be different. For example,
	// for a 8-noded hex, there can be 8 gauss points or there
	// could be 1. The gptr array works like the eptr array, but 
	// allows this difference to occur.
	gptr = (int *)malloc((nelements+1) * sizeof(int));

   int counter = 0;
   int counter0 = 0;
   for (int i = 0; i < nelements; i++) {
	   //printf("(e.%d) - eptr:[%d->%d] - ", i, eptr[i], eptr[i+1]);
	   counter0 = counter;
	   //counter = 0;
	   if (strcmp(ElementType[i], "C3D8") == 0) {
		   //GuassPoints per element
		   GaussPoints[i] = 8;
		   nShapeFunctions[i] = 8;

		   // shp function array needs to hold 8 
		   // shp functions for each of these 8 gauss points
		   // for this one element there are 8 gauss points, 
		   // which each have 8 components in shp function array
		   // so for this element I need 8 * 8 positions to hold
		   counter = counter + nShapeFunctions[i];
	   } 
	   if (strcmp(ElementType[i], "C3D4") == 0) {
		   GaussPoints[i] = 1;
		   nShapeFunctions[i] = 4;
		   // same argument as above
		   counter = counter + nShapeFunctions[i];
	   }
	   //printf("gptr:[%d %d]\n",counter0, counter);
	   gptr[i] = counter0;
	   gptr[i + 1] = counter;
   }
   
   // for debugging purposes
   for (int i = 0; i < nelements; i++) {
	   printf("(e.%d) - eptr:[%d->%d] - [%d->%d]\n", i, eptr[i], eptr[i + 1], gptr[i], gptr[i + 1]);
   }
   printf("size of shp array = %d \n", counter);

   /*set size of shp array  - this holds shp functions for all elements */
   shp =  (double *)malloc(counter*sizeof(double));

   /* initalize shp array with zoros*/
   for (int i = 0; i < nelements; i++) {
	   printf("shp array %d: %d | ", i, nShapeFunctions[i]);
		   for (int k = 0; k < nShapeFunctions[i]; k++) {
	           printf(" %d",gptr[i]+k);
			   shp[gptr[i] + k] = 0.0;
		   }
	   printf("\n");
   }

   for (int i = 0; i < nelements; i++) {
		// Depending on element type call correct shape function library
		// 3D 8-noded hex shape function routine
		if (strcmp(ElementType[i], "C3D8") == 0) {
			double *Chi = (double*)malloc(GaussPoints[i] * ndim * sizeof(double));
			double *GaussWeights = (double*)malloc(GaussPoints[i] * sizeof(double));
			GaussQuadrature3D(i, GaussPoints[i], Chi, GaussWeights);
			for (int k = 0; k < GaussPoints[i]; k++) {
				ShapeFunction_C3D8(i, k, GaussPoints[i], Chi);
			}
			free(Chi);
			free(GaussWeights);
		}
#if 1
		// 3D 4-noded tet shape function routine
		if (strcmp(ElementType[i], "C3D4") == 0) {
			double *Chi = (double*)malloc(GaussPoints[i]* ndim * sizeof(double));
			double *GaussWeights = (double*)malloc(GaussPoints[i] * sizeof(double));
			GaussQuadrature3D(i, GaussPoints[i], Chi, GaussWeights);
			//printf("chi, eta, iota = %f, %f, %f\n", Chi[0], Chi[1], Chi[2]);
			for (int k = 0; k < GaussPoints[i]; k++) {
				ShapeFunction_C3D4(i, k, GaussPoints[i], Chi);
			}
			free(Chi);
			free(GaussWeights);
		}	  
#endif 
   }// loop on nelements
   
   for (int i = 0; i < nelements; i++) {
	   printf("shp array %d: %d | ", i, GaussPoints[i]);
	   for (int k = 0; k < nShapeFunctions[i]; k++) {
		   printf(" %.3f", shp[gptr[i] + k]);
		   //shp[gptr[i] + k] = 0.0;
	   }
	   printf("\n");
   }

   return;
}
