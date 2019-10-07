#include "FemTech.h"

void StressUpdate(int e, int gp){
  FILE_LOG_SINGLE(DEBUGLOGIGNORE, "---- Calling Stress Update %d, %d ---", e, gp);
	//Compressible Neohookean
	if(materialID[pid[e]]==1){
		CompressibleNeoHookean(e, gp);
	}
	// St. Venant-Kirchhoff
	if(materialID[pid[e]]==2){
    FILE_LOG_SINGLE(DEBUGLOGIGNORE, "---- Calling Material Model St. Venant %d, %d ---", e, gp);
		StVenantKirchhoff(e, gp);
	}
  // Linear Elastic
  if(materialID[pid[e]] == 3) {
    LinearElastic(e, gp);
  }
  return;
}
