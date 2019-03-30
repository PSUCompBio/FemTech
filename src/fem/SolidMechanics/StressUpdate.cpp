#include "FemTech.h"

void StressUpdate(int e, int gp){

	//Compressible Neohookean
	if(materialID[e]==1){
		CompressibleNeoHookean(e,gp);
	}
	// St. Venant-Kirchhoff
	if(materialID[e]==2){
		StVenantKirchhoff(e,gp);
	}

	return;
}
