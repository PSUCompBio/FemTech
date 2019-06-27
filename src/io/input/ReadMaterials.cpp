#include "FemTech.h"

void ReadMaterials() {
  /* Function to be called only after reading the mesh */

  FILE *File;
  if ((File = fopen("materials.dat", "r")) == NULL) {
    printf("\nERROR( proc %d ): Cannot open materials.dat file.\n", world_rank);
    exit(0);
  } else {
    int partID;
    for (int i = 0; i < nPIDglobal; i++) {
      fscanf(File, "%d", &partID);
      // printf("elt ID %d\n", elementID);
      fscanf(File, "%d", &materialID[partID]);
      // printf("mat ID %d\n", materialID[elementID]);
      int index = partID * MAXMATPARAMS;
      if (materialID[partID] == 1) { // CompressibleNeoHookean
        // properties[0] = density
        // properties[1] = mu
        // properties[2] = lambda
        fscanf(File, "%lf %lf %lf", &properties[index + 0],
               &properties[index + 1], &properties[index + 2]);
        // printf("element %d Compressible NeoHookean properties (rho, mu, "
        //        "lambda) = %3.3f %3.3f %3.3f\n",
               // elementID, properties[index + 0], properties[index + 1],
               // properties[index + 2]);
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
      } else {
        printf("Material ID for Part %d NOT FOUND! Bye.\n", partID);
        exit(0);
      }
    }
    fclose(File);
  } // end else
  return;
}
