#include "FemTech.h"

void free1DArray(void *array) {
  if(array) {
    free(array);
    array = NULL;
  }
}
// void free2DArray(void *array, int count) {
//   for (int i = 0; i < count; ++i) {
//     if (array[i]) {
//       free(array[i]);
//     }
//   }
//   free1DArray(array);
// }

void FreeArrays() {
  free1DArray(coordinates);
  free1DArray(connectivity);
  free1DArray(pid);
  free1DArray(eptr);
  free1DArray(shp);
  free1DArray(dshp);
  free1DArray(dsptr);
  free1DArray(gptr);
  free1DArray(nShapeFunctions);
  free1DArray(C);
  free1DArray(gaussWeights);
  free1DArray(gpPtr);
  free1DArray(detJacobian);
  free1DArray(mass);
  free1DArray(stiffness);
  free1DArray(nodeDOFDispBC);
  free1DArray(nodeValueDispBC);
  free1DArray(rhs);
	free1DArray(displacements);
  free1DArray(displacement);
  if (ElementType != NULL){
    for (int i = 0; i < nelements; i++){
        free(ElementType[i]);
    }
    free(ElementType);
    ElementType = NULL;
  }
}
