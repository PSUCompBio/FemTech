#include "FemTech.h"

void ReadMaterials() {
  /* Function to be called only after reading the mesh */
  FILE *File;
  if ((File = fopen("materials.dat", "r")) == NULL) {
    printf("\nERROR( proc %d ): Cannot open materials.dat file.\n", world_rank);
    exit(0);
  } 
  int partID;
  printf("DEBUG(%d) : Number of parts = %d\n", world_rank, nPIDglobal);
  bool *checkFullRead = (bool*)malloc(nPIDglobal*sizeof(bool));
  if (checkFullRead == NULL) {
    printf("ERROR(%d) : Error allocating checkFullRead\n", world_rank);
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
      // printf("Part %d Compressible NeoHookean properties (rho, mu, "
      //        "lambda) = %3.3f %3.3f %3.3f\n",
      //         partID, properties[index + 0], properties[index + 1],
      //         properties[index + 2]);
    } else if (materialID[partID] == 2) { // St. Venant
      // properties[0] = density
      // properties[1] = mu
      // properties[2] = lambda
      fscanf(File, "%lf %lf %lf", &properties[index + 0],
              &properties[index + 1], &properties[index + 2]);
      // printf("element %d St. Venant-Kirchhoff properties (rho, mu, lambda) = "
      //        "%3.3f %3.3f %3.3f\n",
              // elementID, properties[index + 0], properties[index + 1],
              // properties[index + 2]);
    } else if (materialID[partID] == 3) { // Linear Elastic
      // properties[0] = density
      // properties[1] = mu
      // properties[2] = lambda
      fscanf(File, "%lf %lf %lf", &properties[index + 0],
              &properties[index + 1], &properties[index + 2]);
      // printf("element %d Linear Elastic properties (rho, mu, lambda) = "
      //        "%3.3f %3.3f %3.3f\n",
              // elementID, properties[index + 0], properties[index + 1],
              // properties[index + 2]);
    } else if (materialID[partID] == 4) { // HGO Isotropic model
      // properties[0] = density
      // properties[1] = mu
      // properties[2] = lambda
      // properties[4] = k1
      // properties[5] = k2
      fscanf(File, "%lf %lf %lf %lf %lf", &properties[index + 0],
              &properties[index + 1], &properties[index + 2], &properties[index + 3], &properties[index + 4]);
      // printf("Part %d HGO properties (rho, mu, lambda, k1, k2) = "
      //        "%3.3f %3.3f %3.3f %3.3f %3.3f\n",
      //         partID, properties[index + 0], properties[index + 1],
      //         properties[index + 2], properties[index + 3], properties[index + 4]);
    } else {
      printf("ERROR : Material ID for Part %d NOT FOUND!\n", partID);
      exit(0);
    }
    checkFullRead[partID] = true;
  }
  fclose(File);
  // Check if all parts have corresponding material properties
  for (int i = 0; i < nPIDglobal; ++i) {
    if(checkFullRead[i] == false) {
      printf("ERROR(%d) : Material properties of partID %d missing\n", world_rank, i);
      exit(0);
    }
  }
  free(checkFullRead);
  
  return;
}
