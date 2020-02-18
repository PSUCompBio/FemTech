#include "utilities.h"

void free1DArray(void *array) {
  if(array) {
    free(array);
    array = NULL;
  }
}
