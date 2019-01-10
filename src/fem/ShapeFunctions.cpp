#include "digitalbrain.h"
void ShapeFunctions(){
   printf("Rank %d: Hello Shape Function!!\n", world_rank);
   printf("\nConnectivity array in processor %d after partitioning =", world_rank);

   for (int i = 0; i < nelements; i++) {
	   printf(" (%d)  ", i);
	   for (int j = eptr[i]; j < eptr[i + 1]; j++) {
		   printf("%d ", connectivity[j]);
	   }
   }
   printf("\n");


   for (int i = 0; i < nelements; i++) {
	   //printf(" (%d)  ", i);
	   //for (int j = eptr[i]; j < eptr[i + 1]; j++) {
		//   printf("%d ", connectivity[j]);
	  // }

	   // Depending on element type call correct shape function library

	   // Call 3D 8-noded hex shape function routine
	   if (strcmp(ElementType[i], "C3D8") == 0) {
		   printf("Computing C3D8 Shape Function , element %d...\n", i);
			 int QuadratureRule = 2;
			 double sf = GaussQuadrature3D(QuadratureRule); /* not sure if this is where it goes*/

			 ShapeFunction_C3D8(sf);

	   }
   }
   printf("\n");


   return;
}
