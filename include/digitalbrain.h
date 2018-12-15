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

// Describes single node that contains ID and coordinates
typedef struct
{
    int Id;
    double XYZ[3];
} NODE;

// Element that contains "Count" nodes
typedef struct
{
    int Id;
    int PId;
    int *Nodes;
    int Count;
    int PartNumber; 
} ELEMENT;

// Partitioning output. Size of "part" is "part_size" 
typedef struct
{
    idx_t edgecut;
    idx_t *part;
    int part_size;
} PARTOUTPUT;

// Partition that contains "Count" elements
typedef struct
{
    ELEMENT *Elements;
    int Count;
    idx_t *Eptr;
} PARTITION;

//-------------------------------------------------------------------------------------------
void MPI_Initialize();

int LineToArray(
    const bool IntOrFloat,
    const bool CheckLastVal,
    const int ColumnToStart,
    const int ColumnCount,
    const char *ConstLine,
    void **Array);

bool ReadInputFile(
    const char *FileName,
    idx_t **elmdist,
    idx_t **MyEptr,
    idx_t **MyEind,
    int   *partArraySize);

bool PartitionMesh(
    idx_t *elmdist,
    idx_t *MyEptr,
    idx_t *MyEind,
    const idx_t NParts,
    const int partArraySize,
    idx_t *edgecut,
    idx_t **part);

PARTOUTPUT *SendReceivePartitioningOutput(
    const idx_t edgecut,
    idx_t *part,
    const int partArraySize);

int CreatePartitions(
    const char *FileName,
    int NParts,
    PARTOUTPUT *Parts,
    PARTITION **PartitionsP);

bool PrepareVTUData(
    const char *FileName,
    const int PartIdx,
    const int NParts,
    const PARTITION Partition);

void WriteVTU(
    const char* FileName,
    const int PartIdx);

void FreeArrays();
void updateConnectivityGlobalToLocal(void);
