#include "FemTech.h"
#include "utilities.h"

void FreeArrays() {
  free1DArray(coordinates);
  free1DArray(connectivity);
  free1DArray(globalNodeID);
  free1DArray(pid);
  free1DArray(global_eid);
  free1DArray(eptr);
  free1DArray(shp);
  free1DArray(dshp);
  free1DArray(dsptr);
  free1DArray(gptr);
  free1DArray(nShapeFunctions);
  free1DArray(C);
  free1DArray(gaussWeights);
  free1DArray(gpPtr);
  free1DArray(detJacobian);
  free1DArray(mass);
  free1DArray(stiffness);
  free1DArray(rhs);
	free1DArray(displacements);
	free1DArray(velocities);
	free1DArray(accelerations);
	free1DArray(accelerations_prev);
	free1DArray(boundary);
	free1DArray(velocities_half);
  free1DArray(fe);
  free1DArray(fi);
	free1DArray(f_net);
	free1DArray(fr_prev);
	free1DArray(fr_curr);
	free1DArray(fi_prev);
	free1DArray(f_damp_prev);
	free1DArray(f_damp_curr);
	free1DArray(displacements_prev);
  free1DArray(F);
	free1DArray(pk2);
	free1DArray(pk2ptr);
	free1DArray(materialID);
	free1DArray(properties);
	free1DArray(detF);
	free1DArray(invF);
	free1DArray(Eavg);
	free1DArray(detFptr);
	free1DArray(InternalsPtr);
	free1DArray(internals);
  // Free arrays used for communication
  free1DArray(recvNodeDisplacement);

  free1DArray(sendProcessID);
  free1DArray(sendNeighbourCount);
  free1DArray(sendNeighbourCountCum);
  free1DArray(sendNodeIndex);
  free1DArray(sendNodeDisplacement);
  free1DArray(stepTime);

  // Temporary arrays used in materials folder
  free1DArray(mat1);
  free1DArray(mat2);
  free1DArray(mat3);
  free1DArray(mat4);
  // InternalForceUpdate
  free1DArray(fintGQ);
  free1DArray(B);
  // Viscoelastic materials
  for (unsigned int i = 0; i < nMax; ++i) {
    free1DArray(Hn[i]);
  }
  free1DArray(Hn);
  free1DArray(S0n);

  if (ElementType != NULL) {
    for (int i = 0; i < nelements; i++){
        free(ElementType[i]);
    }
    free(ElementType);
    ElementType = NULL;
  }

  if (world_rank == 0) {
    fclose(energyFile);
  }
}
