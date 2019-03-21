#include "FemTech.h"
#include "blas.h"

void StressUpdate(int e, int gp){

	// just set material ID to 1 for now  - needs to be removed
	materialID[e]=1;

	//Compressible Neohookean
	if(materialID[e]==1){
		CompressibleNeoHookean(e,gp);
	}

	return;
}
