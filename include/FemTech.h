#ifndef FEMTECH_H
#define FEMTECH_H

#include <stdio.h>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <cstring>
#include <cstdlib>
#include <math.h>
#include <limits.h>

//#include <thread>
//#include <mutex>
//#include <cmath>

#include "GlobalVariables.h"

#include "mpi.h"
#include "parmetis.h"

//-------------------------------------------------------------------------------------------
#define MAX_FILE_LINE 260
#define MAX_ELEMENT_TYPE_SIZE  10

int LineToArray(const bool IntOrFloat, const bool CheckLastVal, \
    const int ColumnToStart, const int ColumnCount, const char *ConstLine, \
    const char *Delim = " \t", void **Array = NULL);
bool ReadInputFile(const char *FileName);
bool PartitionMesh();
void GaussQuadrature3D(int element, int nGaussPoint, double *Chi,double *GaussWeights);
void ShapeFunctions();
void ShapeFunction_C3D8(int e, int gp, double *Chi, double *detJ);
void ShapeFunction_C3D4(int e, int gp, double *Chi, double *detJ);
void ShapeFunction_T3D2(int e, int gp, double *Chi, double *detJ);
void CreateLinearElasticityCMatrix();
void ReadBoundaryCondition(void);
void AllocateArrays();

void Assembly(char *operation);
void StiffnessElementMatrix(double* Ke, int e);
void MassElementMatrix(double* Me, int e);
void WriteVTU(const char* FileName, int step, double time);
void WritePVD(const char* FileName, int step, double time);
void FreeArrays();
void ReadMaterialProperties();
void ApplySteadyBoundaryConditions(void);
void SolveSteadyImplicit(void);
void SolveUnsteadyNewmarkImplicit(double beta, double gamma, double dt, \
    double timeFinal, char* name);
void LumpMassMatrix(void);
void AssembleLumpedMass(void);
void SolveUnsteadyExplicit(double timeFinal, char* name);
void ExplicitDynamics(double timeFinal, char* name);

void GetForce();
void GetForce_3D();
double StableTimeStep();
double CalculateTimeStep(int e);
void CalculateAccelerations();

void CalculateFR();

void CalculateStrain(int e);
void CalculateDeformationGradient(int e, int gp);
void SumOfDeformationGradient(int e, int gp);
void StrainDisplacementMatrix(int e, int gp, int nI, double *B);
void CompressibleNeoHookean(int e, int gp);
void StressUpdate(int e, int gp);
void DeterminateF(int e, int gp);
void InverseF(int e, int gp);
void InternalForceUpdate(int e, int gp, double *force);
void TrussStressForceUpdate(int e, int gp, double *force);
void ReadMaterials();
void StVenantKirchhoff(int e, int gp);
void LinearElastic(int e, int gp);

void inverse3x3Matrix(double* mat, double* invMat, double* det);
//void MultiplyMatrices(double* a, double* b, int sizeM, double* result);
double tripleProduct(double* s, double* a, double* b);
double normOfCrossProduct(double *a, double *b);
void crossProduct(double* a, double* b, double* result);
double norm3D(double *a);
double dotProduct3D(double *a, double *b);


void CheckEnergy(double time);

/* Functions to calculate characteristic lengths */
double CalculateCharacteristicLength(int e);
double CalculateCharacteristicLength_C3D4(int e);
double CalculateCharacteristicLength_C3D8(int e);

/* Functions to calculate center of mass */
void GetBodyCenterofMass(double *cm);
double CalculateCentroidAndVolume(int e, double *cm);
double CalculateCentroidAndVolume_C3D8(int e, double *cm);
void CalculateCentroid_C3D4(int e, double *cm);

void updateMassMatrixNeighbour(void);

/* Geometry related calculations in math folder*/
double volumeHexahedron(double *coordinates);
double areaHexahedronFace(double *coordinates, const int * const index);
double volumeTetrahedron(double *coordinates);

int compare(const void *a, const void *b);
int unique(int *arr, int n);
#endif
