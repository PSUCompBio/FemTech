#include "FemTech.h"

#include <stdlib.h>

void StressUpdate(int e, int gp){
#ifdef DEBUG
	if(debug && 1==0){
  	printf("---- Calling Stress Update %d, %d ---\n", e, gp);
	}
#endif //DEBUG
  switch (materialID[pid[e]]) {
    case 1 : //Compressible Neohookean
             CompressibleNeoHookean(e, gp);
             break;
    case 2 : // St. Venant-Kirchhoff
             StVenantKirchhoff(e, gp);
             break;
    case 3 : // Linear Elastic
             LinearElastic(e, gp);
             break;
    case 4 : // HGO with isotropic fiber distribution
             HGOIsotropic(e, gp);
             break;
    default : printf("Unknown material type\n");
              exit(EXIT_FAILURE);
  }
  return;
}
