#include "FemTech.h"

void StressUpdate(int e, int gp){
  printf("---- Calling Stress Update %d, %d ---\n", e, gp);
	//Compressible Neohookean
	if(materialID[e]==1){
		CompressibleNeoHookean(e,gp);
	}
	// St. Venant-Kirchhoff
	if(materialID[e]==2){
    printf("---- Calling Material Model St. Venant %d, %d ---\n", e, gp);
    printf("---- ----\n");
		StVenantKirchhoff(e,gp);
	}
  // Linear Elastic
  if(materialID[e] == 3) {
    LinearElastic(e, gp);
  }

	return;
}
