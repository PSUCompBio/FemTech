#ifndef DIGITALBRAIN_H
#define DIGITALBRAIN_H

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
//#include "Mesh.h"
//#include "Materials.h"
//#include "BC.h"
//#include "Constraint.h"
#include "mpi.h"
#include "parmetis.h"

//-------------------------------------------------------------------------------------------
#define MAX_FILE_LINE 260
#define MAX_ELEMENT_TYPE_SIZE  10

int LineToArray(
    const bool IntOrFloat,
    const bool CheckLastVal,
    const int ColumnToStart,
    const int ColumnCount,
    const char *ConstLine,
    const char *Delim = " \t",
    void **Array = NULL);
bool ReadInputFile(const char *FileName);
bool PartitionMesh();
void ShapeFunctions();
void ShapeFunction_C3D8();
void GaussQuadrature3D(int QuadratureRule, double *s);
void WriteVTU(const char* FileName);
void FreeArrays();

#endif
