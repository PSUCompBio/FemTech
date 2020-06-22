#include "FemTech.h"

#include <assert.h>

int *materialID;
double *properties;

void ReadMaterials() {
  /* Function to be called only after reading the mesh */
  FILE *File;
  if ((File = fopen("materials.dat", "r")) == NULL) {
    FILE_LOG_SINGLE(ERROR, "Cannot open materials.dat file");
    TerminateFemTech(3);
  } 
  int partID;
  FILE_LOG(INFO, "Number of parts = %d", nPIDglobal);
  bool *checkFullRead = (bool*)malloc(nPIDglobal*sizeof(bool));
  if (checkFullRead == NULL) {
    FILE_LOG_SINGLE(ERROR, "Error allocating checkFullRead");
    TerminateFemTech(12);
  }
  for (int i = 0; i < nPIDglobal; ++i) {
    checkFullRead[i] = false;
  }
	materialID=(int*)calloc(nPIDglobal,sizeof(int));
  if (!materialID) {
    FILE_LOG_SINGLE(ERROR, "Error in allocating materialID array");
    TerminateFemTech(12);
  }
	properties=(double*)calloc(nPIDglobal*MAXMATPARAMS,sizeof(double));
  if (!properties) {
    FILE_LOG_SINGLE(ERROR, "Error in allocating properties array");
    TerminateFemTech(12);
  }

  for (int i = 0; i < nPIDglobal; i++) {
    int count;
    count = fscanf(File, "%d", &partID);
    assert(count == 1);
    count = fscanf(File, "%d", &materialID[partID]);
    assert(count == 1);
    int index = partID * MAXMATPARAMS;
    switch (materialID[partID]) {
      case 0 :// Rigid body
              // properties[0] = density
              count = fscanf(File, "%lf", &properties[index + 0]);
              assert(count == 1);
              break;
      case 1 ://Compressible Neohookean
              // properties[0] = density
              // properties[1] = mu
              // properties[2] = lambda
              count = fscanf(File, "%lf %lf %lf", &properties[index + 0],
                      &properties[index + 1], &properties[index + 2]);
              assert(count == 3);
              // FILE_LOG_SINGLE(DEBUGLOG, "Part %d Compressible NeoHookean properties (rho, mu, "
              //        "lambda) = %3.3f %3.3f %3.3f\n",
              //         partID, properties[index + 0], properties[index + 1],
              //         properties[index + 2]);
              break;
      case 2 :// St. Venant-Kirchhoff
              // properties[0] = density
              // properties[1] = mu
              // properties[2] = lambda
              count = fscanf(File, "%lf %lf %lf", &properties[index + 0],
                      &properties[index + 1], &properties[index + 2]);
              assert(count == 3);
              // FILE_LOG_SINGLE(DEBUGLOG, "element %d St. Venant-Kirchhoff properties (rho, mu, lambda) = "
              //        "%3.3f %3.3f %3.3f\n",
                      // elementID, properties[index + 0], properties[index + 1],
                      // properties[index + 2]);
              break;
      case 3 :// Linear Elastic
              // properties[0] = density
              // properties[1] = mu
              // properties[2] = lambda
              count = fscanf(File, "%lf %lf %lf", &properties[index + 0],
                      &properties[index + 1], &properties[index + 2]);
              assert(count == 3);
              // FILE_LOG_SINGLE(DEBUGLOG, "element %d Linear Elastic properties (rho, mu, lambda) = "
              //        "%3.3f %3.3f %3.3f\n",
                      // elementID, properties[index + 0], properties[index + 1],
                      // properties[index + 2]);
              break;
      case 4 :// HGO with isotropic fiber distribution
              // properties[0] = density
              // properties[1] = mu
              // properties[2] = lambda
              // properties[3] = k1
              // properties[4] = k2
              count = fscanf(File, "%lf %lf %lf %lf %lf", &properties[index + 0],
                      &properties[index + 1], &properties[index + 2], 
                      &properties[index + 3], &properties[index + 4]);
              assert(count == 5);
              // FILE_LOG_SINGLE(DEBUGLOG, "Part %d HGO properties (rho, mu, lambda, k1, k2) = "
              //        "%3.3f %3.3f %3.3f %3.3f %3.3f\n",
              //         partID, properties[index + 0], properties[index + 1],
              //         properties[index + 2], properties[index + 3], properties[index + 4]);
              break;
      case 5 :// HGO with isotropic fiber distribution with 2 term Prony series for viscoelastic properties
              // properties[0] = density
              // properties[1] = mu
              // properties[2] = lambda
              // properties[3] = k1
              // properties[4] = k2
              // properties[5] = g1
              // properties[6] = t1
              // properties[7] = g2
              // properties[8] = t2
              count = fscanf(File, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &properties[index + 0],
                      &properties[index + 1], &properties[index + 2], 
                      &properties[index + 3], &properties[index + 4],
                      &properties[index + 5], &properties[index + 6],
                      &properties[index + 7], &properties[index + 8]);
              assert(count == 9);
              // FILE_LOG_SINGLE(DEBUGLOG, "Part %d HGO properties (rho, mu, lambda, k1, k2) = "
              //        "%3.3f %3.3f %3.3f %3.3f %3.3f\n",
              //         partID, properties[index + 0], properties[index + 1],
              //         properties[index + 2], properties[index + 3], properties[index + 4],
              //         properties[index + 5], properties[index + 6], properties[index + 7],
              //         properties[index + 8]);
              break;
      case 6 :// Viscoelastic material
              // properties[0] = density
              // properties[1] = K
              // properties[2] = G_0
              // properties[3] = G_\infity
              // properties[4] = \beta
              count = fscanf(File, "%lf %lf %lf %lf %lf", &properties[index + 0],
                      &properties[index + 1], &properties[index + 2], 
                      &properties[index + 3], &properties[index + 4]);
              assert(count == 5);
              // FILE_LOG_SINGLE(DEBUGLOG, "Part %d Viscoelastic properties (rho, K, G_0, G_infinity, beta) = "
              //        "%3.3f %3.3f %3.3f %3.3f %3.3f\n",
              //         partID, properties[index + 0], properties[index + 1],
              //         properties[index + 2], properties[index + 3], properties[index + 4]);
              break;
      default :FILE_LOG_SINGLE(ERROR, "Material ID for Part %d not found", partID);
                TerminateFemTech(3);
    }
    checkFullRead[partID] = true;
  }
  fclose(File);
  // Check if all parts have corresponding material properties
  for (int i = 0; i < nPIDglobal; ++i) {
    if(checkFullRead[i] == false) {
      FILE_LOG_SINGLE(ERROR, "Material properties of partID %d missing", i);
      TerminateFemTech(3);
    }
  }
  free(checkFullRead);
  return;
}
