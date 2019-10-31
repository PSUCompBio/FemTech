#include "FemTech.h"

void StressUpdate(int e, int gp){
#ifdef DEBUG
	if(debug && 1==0){
  	printf("---- Calling Stress Update %d, %d ---\n", e, gp);
	}
#endif //DEBUG
	//Compressible Neohookean
	if(materialID[pid[e]]==1){
		CompressibleNeoHookean(e, gp);
	}
	// St. Venant-Kirchhoff
	if(materialID[pid[e]]==2){
#ifdef DEBUG
		if(debug && 1==0){
    	printf("---- Calling Material Model St. Venant %d, %d ---\n", e, gp);
    	printf("---- ----\n");
		}
#endif //DEBUG
		StVenantKirchhoff(e, gp);
	}
  // Linear Elastic
  if(materialID[pid[e]] == 3) {
    LinearElastic(e, gp);
  }
  return;
}
