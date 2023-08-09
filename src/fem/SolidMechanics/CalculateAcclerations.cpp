#include "FemTech.h"
#include "blas.h"

void CalculateAccelerations(){
  FILE_LOGArray(DEBUGLOGIGNORE, f_net, nDOF, "RHS of Solution");
  if(embed){
    for (int i = 0; i < nDOF; i++) {
	   int embedid = int(i/3);
		if(nodeconstrain[embedid]!=-1)
	      	   accelerations[i] = 0;
		else accelerations[i] = (f_net[i]-dynamicDamping*velocities_half[i])/mass[i];
    }
  }
  else
  for (int i = 0; i < nDOF; ++i) {
    accelerations[i] = (f_net[i]-dynamicDamping*velocities_half[i])/mass[i];
  }
	return;
}
