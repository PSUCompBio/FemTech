#include "FemTech.h"

double CalculateCharacteristicLength(int e) {
  double cl;
  if (strcmp(ElementType[e], "C3D8") == 0) {
    cl = CalculateCharacteristicLength_C3D8(e);
  } else {
    if (strcmp(ElementType[e], "C3D4") == 0) {
      cl = CalculateCharacteristicLength_C3D4(e);
    } else {
	if (strcmp(ElementType[e], "C3D8R") == 0) {
      cl = CalculateCharacteristicLength_C3D8(e);
    } 	else {
      FILE_LOG_SINGLE(ERROR, "Unknown Element Type Encountered");
      cl = 1.0;
    }
  }
}
  return cl;
}
