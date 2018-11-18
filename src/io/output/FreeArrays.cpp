#include "digitalbrain.h"

void FreeArrays(){

free(coordinates);
free(connectivity);
free(mid);
free(pid);
free(eptr);

return;
}
