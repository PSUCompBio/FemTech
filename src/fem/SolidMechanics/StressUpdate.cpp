#include "FemTech.h"

void StressUpdate(int e, int gp){
#ifdef DEBUG
	if(debug){
  	printf("---- Calling Stress Update %d, %d ---\n", e, gp);
	}
#endif //DEBUG
	//Compressible Neohookean
	if(materialID[e]==1){
		CompressibleNeoHookean(e,gp);
	}
	// St. Venant-Kirchhoff
	if(materialID[e]==2){
#ifdef DEBUG
		if(debug){
    	printf("---- Calling Material Model St. Venant %d, %d ---\n", e, gp);
    	printf("---- ----\n");
		}
#endif //DEBUG
		StVenantKirchhoff(e,gp);
	}
  // Linear Elastic
  if(materialID[e] == 3) {
    LinearElastic(e, gp);
  }
  return;
}
