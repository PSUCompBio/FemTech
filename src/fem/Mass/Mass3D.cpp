#include "FemTech.h"

void Mass3D(){
	// set the debug flag for this file
	int debug = 1;

	// for debugging
	if (debug) {
		for (int i = 0; i < nelements; i++) {
			printf("shp array e.%d with %d shp functions\n", i, GaussPoints[i]);
			for (int k = 0; k < nShapeFunctions[i]; k++) {
				printf("int %d | dshp: %8.5f %8.5f %8.5f shp: %8.5f\n", k, dshp[dsptr[i] + k * ndim + 0],
														   dshp[dsptr[i] + k * ndim + 1],
														   dshp[dsptr[i] + k * ndim + 2],
														   shp[gptr[i] + k]);

			}
			printf("\n");
		}
	}
	return;

	return;
}
