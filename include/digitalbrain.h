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

//#include <thread>
//#include <mutex>
//#include <cmath>

#include "GlobalVariables.h"
//#include "Mesh.h"
//#include "Materials.h"
//#include "BC.h"
//#include "Constraint.h"
#if PARALLEL
	#include "mpi.h"
	#include "parmetis.h"
#endif

void MPI_Initialize();

bool ReadInputFile(const char* str);

void PartitionMesh();

void WriteVTU(char* str);

void FreeArrays();
