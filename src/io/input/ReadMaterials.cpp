#include "FemTech.h"

#include <assert.h>

int *materialID;
double *properties;
double *waveSpeed;
bool *viscoElasticPart;

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
  viscoElasticPart=(bool*)malloc(nPIDglobal*sizeof(bool));

  for (int i = 0; i < nPIDglobal; i++) {
    int count;
    count = fscanf(File, "%d", &partID);
    assert(count == 1);
    count = fscanf(File, "%d", &materialID[partID]);
    assert(count == 1);
    int index = partID * MAXMATPARAMS;
    int nTerms = 0;
    int nPronyL = 0;
    switch (materialID[partID]) {
      case 0 :// Rigid body
              // properties[0] = density
              count = fscanf(File, "%lf", &properties[index + 0]);
              assert(count == 1);
              viscoElasticPart[partID] = false;
              break;
      case 1 ://Compressible Neohookean
              // properties[0] = density
              // properties[1] = mu
              // properties[2] = lambda
              count = fscanf(File, "%lf %lf %lf", &properties[index + 0],
                      &properties[index + 1], &properties[index + 2]);
              assert(count == 3);
              viscoElasticPart[partID] = false;
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
              viscoElasticPart[partID] = false;
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
              viscoElasticPart[partID] = false;
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
              viscoElasticPart[partID] = false;
              assert(count == 5);
              // FILE_LOG_SINGLE(DEBUGLOG, "Part %d HGO properties (rho, mu, lambda, k1, k2) = "
              //        "%3.3f %3.3f %3.3f %3.3f %3.3f\n",
              //         partID, properties[index + 0], properties[index + 1],
              //         properties[index + 2], properties[index + 3], properties[index + 4]);
              break;
      case 5 :// HGO with isotropic fiber distribution with 2 term Prony series for viscoelastic properties
              // properties[0] = density
              // properties[1] = mu
              // properties[2] = K
              // properties[3] = k1
              // properties[4] = k2
              // properties[5] = n
              // properties[6+2*i] = gi
              // properties[7+2*i] = ti
              count = fscanf(File, "%lf %lf %lf %lf %lf %lf", &properties[index + 0],
                      &properties[index + 1], &properties[index + 2], 
                      &properties[index + 3], &properties[index + 4],
                      &properties[index + 5]);
              assert(count == 6);
              if (int(properties[index+5]) > 6) {
                FILE_LOG_SINGLE(ERROR, "Material ID 5, n = %d > 6", int(properties[index+5]));
                TerminateFemTech(3);
              }
              for (int j = 0; j < properties[index+5]; ++j) {
                count = fscanf(File, "%lf %lf", &properties[index + 6 + 2*j], 
                    &properties[index + 7 + 2*j]);
                assert(count == 2);
              }
              viscoElasticPart[partID] = true;
              // FILE_LOG_SINGLE(WARNING, "Part %d HGO properties (rho, mu, K, k1, k2, n) = "
              //        "%3.3f %3.3f %3.3f %3.3f %3.3f %3.3f\n",
              //         partID, properties[index + 0], properties[index + 1],
              //         properties[index + 2], properties[index + 3], properties[index + 4],
              //         properties[index + 5]);
              break;
      case 6 :// lsDyna KM Equivalent material
              // properties[0] = density
              // properties[1] = K
              // properties[2] = G_0
              // properties[3] = G_\infity
              // properties[4] = \beta 
              count = fscanf(File, "%lf %lf %lf %lf %lf", &properties[index + 0],
                      &properties[index + 1], &properties[index + 2], 
                      &properties[index + 3], &properties[index + 4]);
              assert(count == 5);
              // Convert G_0 to prony series equivalent g_1 : 
              // g_1 = (G_0 - / G_\infity)/G_\infity = (G_0/G_\infity -1 )
              // Store in properties[3]
              properties[index+3] = properties[index+3]/properties[index+2] - 1.0;
              // Convert properties[1], K to \lambda
              // \lambda = K - 2G/3
              properties[index+1] = properties[index+1] - 2.0*properties[index+2]/3.0;
              viscoElasticPart[partID] = true;
              // FILE_LOG_SINGLE(WARNING, "Part %d Viscoelastic properties (rho, K, G_0, G_infinity, beta) = "
              //        "%3.3f %3.3f %3.3f %3.3f %3.3f\n",
              //         partID, properties[index + 0], properties[index + 1],
              //         properties[index + 2], properties[index + 3], properties[index + 4]);
              // Convert \beta to \tau
              // properties[index+4] = 1.0/properties[index+4];
              break;
      case 7 :// Ogden model with max N = 3
              // properties[0] = density
              // properties[1] = K
              // properties[2] = nOgden : number of terms in ogden series (max 3)
              // properties[3+2*i] = \alpha_i
              // properties[4+2*i] = \mu_i
              count = fscanf(File, "%lf %lf %lf", &properties[index + 0],
                      &properties[index + 1], &properties[index + 2]);
              assert(count == 3);
              nTerms = static_cast<int>(properties[index+2]);
              if (nTerms > 3) {
                FILE_LOG_SINGLE(ERROR, "Material ID 7 (Ogden), Number of Terms = %d > 3", nTerms);
                TerminateFemTech(3);
              }
              for (int j = 0; j < nTerms; ++j) {
                count = fscanf(File, "%lf %lf", &properties[index + 3 + 2*j],
                    &properties[index + 4 + 2*j]);
                assert(count == 2);
              }
              viscoElasticPart[partID] = false;
              break;
      case 8 :// Ogden Viscoelastic model with max N = 3, max N Prony = 6
              /* \rho      = properties(0)
              * K         = properties(1)
              * N_{Ogden} = properties(2), maximum value of 3
              * \alpha_i  = properties(3+2*i) i = 0 to N_{Ogden}-1
              * \mu_i     = properties(4+2*i) i = 0 to N_{Ogden}-1
              * N_{Prony} = properties(3+N_{Ogden}*2), maximum value of 6
              * g_j       = properties(4+2*N_{Ogden}+2*j) j = 0 to N_{Prony}-1
              * \tau_j    = properties(5+2*N_{Ogden}+2*j) j = 0 to N_{Prony}-1
              * MAXMATPARAMS = 4+2*(N_{Ogden}+N_{Prony}) = 22
              * */
              count = fscanf(File, "%lf %lf %lf", &properties[index + 0],
                      &properties[index + 1], &properties[index + 2]);
              assert(count == 3);
              nTerms = static_cast<int>(properties[index+2]);
              if (nTerms > 3) {
                FILE_LOG_SINGLE(ERROR, "Material ID 7 (Ogden), Number of Terms = %d > 3", nTerms);
                TerminateFemTech(3);
              }
              for (int j = 0; j < nTerms; ++j) {
                count = fscanf(File, "%lf %lf", &properties[index + 3 + 2*j],
                    &properties[index + 4 + 2*j]);
                assert(count == 2);
              }
              count = fscanf(File, "%lf", &properties[index + 3 + 2*nTerms]);
              assert(count == 1);
              nPronyL = static_cast<int>(properties[index+3+2*nTerms]);
              for (int j = 0; j < nPronyL; ++j) {
                count = fscanf(File, "%lf %lf", &properties[index + 4 + 2*nTerms + 2*j],
                    &properties[index + 5 + 2*nTerms + 2*j]);
                assert(count == 2);
              }
              viscoElasticPart[partID] = true;
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
  // Compute wave speeds and store for all materials
	waveSpeed = (double*)malloc(nPIDglobal*sizeof(double));
  for (int i = 0; i < nPIDglobal; ++i) {
    waveSpeed[i] = CalculateWaveSpeed(i);
  }
  return;
}
