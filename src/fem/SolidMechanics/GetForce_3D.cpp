#include "FemTech.h"
#include "blas.h"


// this might be able to be moved to top-level getforce function
// as long as dimension and 2d/3D is accounted for.
void GetForce_3D(){
	for (int i=0;i<nelements;i++){
		for(int j=0;j<GaussPoints[i];j++){
			CalculateDeformationGradient(i,j);
			DeterminateF(i,j);
			InverseF(i,j);
			StressUpdate(i,j);
			InternalForceUpdate(i,j);
		} //loop on gauss points
	}// loop on i, nelements
	return;
}
