#include "FemTech.h"

#include <assert.h>

/* Global Variables */
int nparts = 0;
int nelements = 0;    // number of elements in processor
int nallelements = 0; // number of all elements in mesh
int nNodes = 0;
int ndim = 0;
int nDOF = 0;
double *coordinates;
int *connectivity;
int *eptr;
int *pid;
int *global_eid;
char **ElementType;
int nPIDglobal = 0;
int nPID = 0;

//-------------------------------------------------------------------------------------------
int StrCmpCI(const char *str1, const char *str2);
bool ReadLsDyna(const char *FileName);
bool ReadAbaqus(const char *FileName);
//-------------------------------------------------------------------------------------------
void ReadInputFile(const char *FileName) {
  // Checking MPI variables validity
  if (world_size < 1 || world_rank < 0 || world_rank >= world_size) {
    FILE_LOG_SINGLE(ERROR, "'World_size' and/or 'world_rank' variable is not valid");
    TerminateFemTech(3);
  }

  // Checking file name validity
  if (FileName == NULL || strlen(FileName) == 0) {
    FILE_LOG_SINGLE(ERROR, "Input file name is empty");
    TerminateFemTech(3);
  }

  // Checking file extension and calling corresponding reader
  const char *Ext = "";
  for (int i = strlen(FileName) - 1; i >= 0; i--) {
    if (FileName[i] == '.') {
      Ext = &FileName[i];
      break;
    }
  }
  int returnValue = 0;
  if (StrCmpCI(Ext, ".k") == 0) {
    returnValue = ReadLsDyna(FileName);
  } else if (StrCmpCI(Ext, ".inp") == 0) {
    returnValue = ReadAbaqus(FileName);
  } else {
    FILE_LOG_SINGLE(ERROR, "Input file type is unknown");
  }
  if (returnValue) {
    // Change pid from 1 based number to zero based numbering
    // Get ans set unique number of partIDs both local and global
    // Find the unique number of pid in the mesh
    int *sortedPID = (int *)malloc(nelements * sizeof(int));
    memcpy(sortedPID, pid, nelements * sizeof(int));
    qsort(sortedPID, nelements, sizeof(int), compare);
    nPID = unique(sortedPID, nelements);
    int globalPIDmax = sortedPID[nPID - 1] + 1;
    MPI_Allreduce(MPI_IN_PLACE, &globalPIDmax, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);
    nPIDglobal = globalPIDmax;
    nDOF = nNodes * ndim;

    FILE_LOGArrayInt(DEBUGLOGIGNORE, sortedPID, nPID, \
        "Number of local pid : %d, global : %d", nPID, nPIDglobal);

    free(sortedPID);
    assert(nPID <= nPIDglobal);
    assert(nPIDglobal > 0);
  } else {
    FILE_LOG_SINGLE(ERROR, "Input mesh file read failed");
    TerminateFemTech(3);
  }
}
//-------------------------------------------------------------------------------------------
/*
   Converts line into array of int or double if line contains
   valid int/double separated by whitespase or tab.
*/
int LineToArray(
    const bool IntOrFloat, // If true, will be converted to int else to double
    const bool
        CheckLastVal, // If true, if next column's value is the same as the last
                      // column's one, next coulmns won't be fetched
    const int ColumnToStart, // One-based index of column to start to fetch
    const int ColumnCount,   // Number of columns to fetch. Will be used if
                             // greater than zero
    const char *ConstLine,   // Line to convert
    const char *Delim,       // Delimeter used in "strtok()" function
    void **Array) {          // Output array, int or double

  int Result = 0;
  char Line[MAX_FILE_LINE];
  strcpy(Line, ConstLine);
  double LastVal = INT_MAX;
  int Column = 0, LastRepeatedValuesCount = 0;
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
    if (CheckLastVal) {
      if (Val <= 0) {
        LastRepeatedValuesCount = 0;
        break;
      } else if (Val == LastVal) {
        LastRepeatedValuesCount = LastRepeatedValuesCount + 1;
      } else {
        LastRepeatedValuesCount = 0;
      }
    }
    LastVal = Val;
    Result = Result + 1;
  }

  if (Result > 0) {
    if (CheckLastVal && LastRepeatedValuesCount > 1) {
      Result = Result - LastRepeatedValuesCount;
      if (Result < 0) {
        Result = 0;
      }
    }
    if (ColumnCount > 0 && Result > ColumnCount) {
      Result = ColumnCount;
    }
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
  } else {
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
    } else {
      DoubleArray[i] = Val;
    }
    i = i + 1;
  }
  return Result;
}
