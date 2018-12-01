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
// Partitioning output
typedef struct
{
    idx_t edgecut;
    idx_t *part;
    size_t part_size;
} PARTOUTPUT;

// Element that contains "Count" nodes
typedef struct
{
    int *Nodes;
    size_t Count;
} ELEMENT;


//-------------------------------------------------------------------------------------------
void MPI_Initialize();

bool ReadInputFile(
    const char *FileName,
    idx_t **elmdist,
    idx_t **MyEptr,
    idx_t **MyEind,
    size_t *partArraySize);

bool PartitionMesh(
    idx_t *elmdist,
    idx_t *MyEptr,
    idx_t *MyEind,
    const idx_t NParts,
    const size_t partArraySize,
    idx_t *edgecut,
    idx_t **part);

PARTOUTPUT *SendReceivePartitioningOutput(
    const idx_t edgecut,
    idx_t *part,
    const size_t partArraySize);

void WriteVTU(char* str);

void FreeArrays();
