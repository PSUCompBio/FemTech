#include "FemTech.h"

void ReadMaterials() {
  /* Function to be called only after reading the mesh */
  FILE *File;
  if ((File = fopen("materials.dat", "r")) == NULL) {
    FILE_LOG_SINGLE(ERROR, "Cannot open materials.dat file");
    exit(0);
  } 
  int partID;
  FILE_LOG(INFO, "Number of parts = %d", nPIDglobal);
  bool *checkFullRead = (bool*)malloc(nPIDglobal*sizeof(bool));
  if (checkFullRead == NULL) {
    FILE_LOG_SINGLE(ERROR, "Error allocating checkFullRead");
    exit(0);
  }
  for (int i = 0; i < nPIDglobal; ++i) {
    checkFullRead[i] = false;
  }

  for (int i = 0; i < nPIDglobal; i++) {
    fscanf(File, "%d", &partID);
    fscanf(File, "%d", &materialID[partID]);
    int index = partID * MAXMATPARAMS;
    if (materialID[partID] == 1) { // CompressibleNeoHookean
      // properties[0] = density
      // properties[1] = mu
      // properties[2] = lambda
      fscanf(File, "%lf %lf %lf", &properties[index + 0],
              &properties[index + 1], &properties[index + 2]);
    } else if (materialID[partID] == 2) { // St. Venant
      // properties[0] = density
      // properties[1] = mu
      // properties[2] = lambda
      fscanf(File, "%lf %lf %lf", &properties[index + 0],
              &properties[index + 1], &properties[index + 2]);
    } else if (materialID[partID] == 3) { // Linear Elastic
      // properties[0] = density
      // properties[1] = mu
      // properties[2] = lambda
      fscanf(File, "%lf %lf %lf", &properties[index + 0],
              &properties[index + 1], &properties[index + 2]);
    } else {
      FILE_LOG_SINGLE(ERROR, "Material ID for Part %d not found", partID);
      exit(0);
    }
    checkFullRead[partID] = true;
  }
  fclose(File);
  // Check if all parts have corresponding material properties
  for (int i = 0; i < nPIDglobal; ++i) {
    if(checkFullRead[i] == false) {
      FILE_LOG_SINGLE(ERROR, "Material properties of partID %d missing", i);
      exit(0);
    }
  }
  free(checkFullRead);
  return;
}
