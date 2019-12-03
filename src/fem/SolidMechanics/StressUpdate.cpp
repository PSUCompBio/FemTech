#include "FemTech.h"

#include <stdlib.h>

void StressUpdate(int e, int gp){
#ifdef DEBUG
	if(debug && 1==0){
  	printf("---- Calling Stress Update %d, %d ---\n", e, gp);
	}
#endif //DEBUG
  switch (materialID[pid[e]]) {
	           //Compressible Neohookean
    case 1 : CompressibleNeoHookean(e, gp);
             break;
	           // St. Venant-Kirchhoff
    case 2 : StVenantKirchhoff(e, gp);
             break;
             // Linear Elastic
    case 3 : LinearElastic(e, gp);
             break;
             // HGO with isotropic fiber distribution
    case 4 : HGOIsotropic(e, gp);
             break;
    default : printf("Unknown material type\n");
              exit(EXIT_FAILURE);
  }
  return;
}
