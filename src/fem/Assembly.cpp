#include "FemTech.h"

double *mass;

void Assembly(char *operation){
	// set the debug flag for this file
	int debug = 0;
	
	if (strcmp("mass", operation) == 0) {
		
		Mass3D();

	}

	return;
}
