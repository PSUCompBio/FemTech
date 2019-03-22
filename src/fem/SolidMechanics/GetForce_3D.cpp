#include "FemTech.h"
#include "blas.h"

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
