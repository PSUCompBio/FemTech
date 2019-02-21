#include "FemTech.h"

/* Global Variables */
int nparts=0;
int nelements=0;        // number of elements in processor
int nallelements=0;     // number of all elements in mesh
int nnodes=0;
int ndim=0;
double *coordinates;
int *connectivity;
int *eptr;
int *pid;
int *mid;
char **ElementType;
int world_rank;
int world_size;

//-------------------------------------------------------------------------------------------
int StrCmpCI(const char *str1, const char *str2);
bool ReadLsDyna(const char *FileName);
bool ReadAbaqus(const char *FileName);
//-------------------------------------------------------------------------------------------
bool ReadInputFile(const char *FileName) {
    // Checking MPI variables validity
    if (world_size < 1 || world_rank < 0 || world_rank >= world_size) {
        printf("\nERROR( proc %d ): 'world_size' and/or 'world_rank' variable is not valid.\n", world_rank);
        return false;
    }

    // Checking file name validity
    if (FileName == NULL || strlen(FileName) == 0) {
        printf("\nERROR( proc %d ): Input file name is empty.\n", world_rank);
        return false;
    }
    
    // Checking file extension and calling corresponding reader
    const char *Ext = "";
    for (int i = strlen(FileName) - 1; i >= 0; i--) {
        if (FileName[i] == '.') {
            Ext = &FileName[i];
            break;
        }
    }
    if (StrCmpCI(Ext, ".k") == 0) {
        return ReadLsDyna(FileName);
    }
    else if (StrCmpCI(Ext, ".inp") == 0) {
        return ReadAbaqus(FileName);
    }
    else {
        printf("\nERROR( proc %d ): Input file type is unknown.\n", world_rank);
    }

    return false;
}
//-------------------------------------------------------------------------------------------
/*
   Converts line into array of int or double if line contains
   valid int/double separated by whitespase or tab.
*/
int LineToArray(
    const bool IntOrFloat,      // If true, will be converted to int else to double
    const bool CheckLastVal,    // If true, if next column's value is the same as the last column's one, next coulmns won't be fetched
    const int ColumnToStart,    // One-based index of column to start to fetch
    const int ColumnCount,      // Number of columns to fetch. Will be used if greater than zero
    const char *ConstLine,      // Line to convert
    const char *Delim,          // Delimeter used in "strtok()" function
    void **Array) {             // Output array, int or double

    int Result = 0;
    char Line[MAX_FILE_LINE];
    strcpy(Line, ConstLine);
    double LastVal = INT_MAX;
    int Column = 0;
    char *T = strtok(Line, Delim);
    while (T != NULL) {
        char *EndPtr;
        const double Val = strtod(T, &EndPtr);
        if (EndPtr == T) {
            // Error. Token in a line cannot be converted to number.
            return 0;
        }
        T = strtok(NULL, Delim);
        Column = Column + 1;
        if (Column < ColumnToStart) {
            // Skip columns to go to required ones.
            continue;
        }
        if ((CheckLastVal && (Val == 0 || Val == LastVal)) || (ColumnCount > 0 && Result >= ColumnCount)) {
            // Stop fetching columns.
            break;
        }
        LastVal = Val;
        Result = Result + 1;
    }

    if (Result == 0 || Array == NULL) {
        return Result;
    }

    strcpy(Line, ConstLine);
    int *IntArray;
    double *DoubleArray;
    if (IntOrFloat) {
        IntArray = (int *)calloc(Result, sizeof(int));
        *Array = (void *)IntArray;
    }
    else {
        DoubleArray = (double *)calloc(Result, sizeof(double));
        *Array = (void *)DoubleArray;
    }
    
    int i = 0;
    Column = 0;
    T = strtok(Line, Delim);
    while (T != NULL && i < Result) {
        const double Val = strtod(T, NULL);
        T = strtok(NULL, Delim);
        Column = Column + 1;
        if (Column < ColumnToStart) {
            continue;
        }
        if (IntOrFloat) {
            IntArray[i] = int(Val);
        }
        else {
            DoubleArray[i] = Val;
        }
        i = i + 1;
    }
    return Result;
}

int *nodeDOFDispBC;
double *nodeValueDispBC;
int nSpecifiedDispBC = 0;

void ReadBoundaryCondition(void) {
  // TODO(Anil) Convert hardcoded BC to that from file inputs
  nSpecifiedDispBC = 10;
  nodeDOFDispBC = (int*)malloc(nSpecifiedDispBC*sizeof(int));
  nodeValueDispBC = (double*)malloc(nSpecifiedDispBC*sizeof(double));
  // Node one x, y, z displacements set to zero
  nodeDOFDispBC[0] = 0;
  nodeValueDispBC[0] = 0;
  nodeDOFDispBC[1] = 1;
  nodeValueDispBC[1] = 0;
  nodeDOFDispBC[2] = 2;
  nodeValueDispBC[2] = 0;

  nodeDOFDispBC[3] = 5;
  nodeValueDispBC[3] = 0;
  nodeDOFDispBC[4] = 8;
  nodeValueDispBC[4] = 0.1;
  nodeDOFDispBC[5] = 11;
  nodeValueDispBC[5] = 0.1;
  nodeDOFDispBC[6] = 14;
  nodeValueDispBC[6] = 0;
  nodeDOFDispBC[7] = 17;
  nodeValueDispBC[7] = 0;
  nodeDOFDispBC[8] = 20;
  nodeValueDispBC[8] = 0.1;
  nodeDOFDispBC[9] = 23;
  nodeValueDispBC[9] = 0.1;
}
