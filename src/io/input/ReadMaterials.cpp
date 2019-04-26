#include "FemTech.h"

void ReadMaterials() {
  FILE *File;
  if ((File = fopen("materials.dat", "r")) == NULL) {
    printf("\nERROR( proc %d ): Cannot open materials.dat file.\n", world_rank);
    exit(0);
  } else {
    int elementID;
    for (int i = 0; i < nelements; i++) {
      fscanf(File, "%d", &elementID);
      // printf("elt ID %d\n", elementID);
      fscanf(File, "%d", &materialID[elementID]);
      // printf("mat ID %d\n", materialID[elementID]);
      int index = elementID * MAXMATPARAMS;
      if (materialID[elementID] == 1) { // CompressibleNeoHookean
        // properties[0] = density
        // properties[1] = mu
        // properties[2] = lambda
        fscanf(File, "%lf %lf %lf", &properties[index + 0],
               &properties[index + 1], &properties[index + 2]);
        // printf("element %d Compressible NeoHookean properties (rho, mu, "
        //        "lambda) = %3.3f %3.3f %3.3f\n",
               // elementID, properties[index + 0], properties[index + 1],
               // properties[index + 2]);
      } else if (materialID[elementID] == 2) { // St. Venant
        // properties[0] = density
        // properties[1] = mu
        // properties[2] = lambda
        fscanf(File, "%lf %lf %lf", &properties[index + 0],
               &properties[index + 1], &properties[index + 2]);
        // printf("element %d St. Venant-Kirchhoff properties (rho, mu, lambda) = "
        //        "%3.3f %3.3f %3.3f\n",
               // elementID, properties[index + 0], properties[index + 1],
               // properties[index + 2]);
      } else if (materialID[elementID] == 3) { // Linear Elastic
        // properties[0] = density
        // properties[1] = mu
        // properties[2] = lambda
        fscanf(File, "%lf %lf %lf", &properties[index + 0],
               &properties[index + 1], &properties[index + 2]);
        // printf("element %d Linear Elastic properties (rho, mu, lambda) = "
        //        "%3.3f %3.3f %3.3f\n",
               // elementID, properties[index + 0], properties[index + 1],
               // properties[index + 2]);
      } else {
        printf("Material ID for Element %d NOT FOUND! Bye.\n", elementID);
        exit(0);
      }
    }
    fclose(File);
  } // end else
  return;
}
