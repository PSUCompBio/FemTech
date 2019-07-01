#include "FemTech.h"
#include "blas.h"

/*! \brief Calculates the resultant nodal force after each time step.
 *
 * This function represents the 'getforce' step in Belytschko (Box 6.1 - Explicit FEM Algorithm).
 *  For each hex element, this function calculates the internal nodal force vector
 * and the resultant nodal force vector. Once, this is calculated for each element,
 * the resultant vectors are scattered into global vectors.
 */

void GetForce() {

	if(ndim == 1){
		GetForce_3D();
		printf("Here we are outside of the gates of valhala\n");
		exit (EXIT_FAILURE);
	}
	else if (ndim == 2){
		printf("GetForce function not yet implemented for 2D. Bye.\n");
		exit(0);
	}
	else if (ndim ==3){
		GetForce_3D();
  }
}
