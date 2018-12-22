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
#define MAX_FILE_LINE 128
#define MAX_ELEMENT_TYPE_SIZE  10

void MPI_Initialize();

int LineToArray(
    const bool IntOrFloat,
    const bool CheckLastVal,
    const int ColumnToStart,
    const int ColumnCount,
    const char *ConstLine,
    void **Array);

bool ReadInputFile(const char *FileName);

bool PartitionMesh();

void WriteVTU(const char* FileName);

void FreeArrays();
//void updateConnectivityGlobalToLocal(void);
