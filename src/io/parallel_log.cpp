#include "parallel_log.h"
#include "GlobalVariables.h"

#include <string.h>
#include <stdio.h>
#include <stdarg.h>

#include "mpi.h"

MPI_File logFilePtr;
// Char buffer for simple text log
char s1[256];
char s2[256];
va_list arg;

const char* levelToString(enum logLevel level) {
	static const char* const buffer[] = {"ERROR  ", "WARNING", "INFO   ", "DEBUG  "};
  return buffer[level];
}

void initLog(const char *outFileName) {
  MPI_Info infoin;
  MPI_Info_create(&infoin);
  MPI_Info_set(infoin, "access_style", "write_once,random");

  int err;
  err = MPI_File_open(MPI_COMM_WORLD, outFileName, MPI_MODE_EXCL|MPI_MODE_WRONLY|MPI_MODE_CREATE, infoin, &logFilePtr);
  if (err != MPI_SUCCESS) {
    if (world_rank == 0) {
      MPI_File_delete(outFileName, MPI_INFO_NULL);
    }
    err = MPI_File_open(MPI_COMM_WORLD, outFileName, MPI_MODE_EXCL|MPI_MODE_WRONLY|MPI_MODE_CREATE, infoin, &logFilePtr);
    if (err != MPI_SUCCESS) {
      fprintf(stdout, "ERROR(%d): Unable to open output log file\n", world_rank);
    }
  }
}

void finaliseLog() {
  MPI_File_close(&logFilePtr);
}

void fileLog(enum logLevel level, const char* fmt, ...) {
  sprintf(s1, "%s(%d): %s\n", levelToString(level), world_rank, fmt);
  va_start(arg, fmt);
  vsprintf(s2, s1, arg);
  va_end(arg);
  MPI_File_write_ordered(logFilePtr, s2, strlen(s2), MPI_CHAR, MPI_STATUS_IGNORE);
}

void fileLogMaster(enum logLevel level, const char* fmt, ...) {
  if (world_rank == 0) {
    sprintf(s1, "%s(%d): %s\n", levelToString(level), world_rank, fmt);
    va_start(arg, fmt);
    vsprintf(s2, s1, arg);
    va_end(arg);
    MPI_File_write_shared(logFilePtr, s2, strlen(s2), MPI_CHAR, MPI_STATUS_IGNORE);
  }
}

void fileLogSingle(enum logLevel level, const char* fmt, ...) {
  sprintf(s1, "%s(%d): %s\n", levelToString(level), world_rank, fmt);
  va_start(arg, fmt);
  vsprintf(s2, s1, arg);
  va_end(arg);
  MPI_File_write_shared(logFilePtr, s2, strlen(s2), MPI_CHAR, MPI_STATUS_IGNORE);
}
// void fileLogMatrix(enum logLevel level, const double* mat, const int n, \
//     const int m, const char* txt) {
//   int world_rank = 1;
//   fprintf(fptr, "%s (%d) : %s\n", levelToString(level), world_rank, txt);
//   for (int i = 0; i < n; ++i) {
//     for (int j = 0; j < m; ++j) {
//       fprintf(fptr, "%12.8e  ", mat[i+j*n]);
//     }
//     fprintf(fptr, "\n");
//   }
// }
//
// void fileLogArray(enum logLevel level, const double* arr, const int n, \
//     const char* txt) {
//   int world_rank = 1;
//   fprintf(fptr, "%s (%d) : %s\n", levelToString(level), world_rank, txt);
//   for (int i = 0; i < n; ++i) {
//     fprintf(fptr, "%12.8e\n", arr[i]);
//   }
// }
//
// void fileLogArrayPartial(enum logLevel level, const double* arr, const int nS, \
//     const int nE, const char* txt) {
//   int world_rank = 1;
//   fprintf(fptr, "%s (%d) : %s\n", levelToString(level), world_rank, txt);
//   for (int i = nS; i < nE; ++i) {
//     fprintf(fptr, "%12.8e\n", arr[i]);
//   }
// }
//
// void fileLogMatrixPartial(enum logLevel level, const double* mat, const int n, \
//     const int nS, const int nE, const int mS, const int mE, const char* txt) {
//   int world_rank = 1;
//   fprintf(fptr, "%s (%d) : %s\n", levelToString(level), world_rank, txt);
//   for (int i = nS; i < nE; ++i) {
//     for (int j = mS; j < mE; ++j) {
//       fprintf(fptr, "%12.8e  ", mat[i+j*n]);
//     }
//     fprintf(fptr, "\n");
//   }
// }


// #define FILE_LOG(level, ...) \
//     if (level > FILELOG_MAX_LEVEL) ;\
//     else fileLog(level, __VA_ARGS__);
//
// #define FILE_LOGMatrix(level, mat, n, m, txt) \
//     if (level > FILELOG_MAX_LEVEL) ;\
//     else fileLogMatrix(level, mat, n, m, txt);
//
// #define FILE_LOGMatrixPartial(level, mat, n, nS, nE, mS, mE, txt) \
//     if (level > FILELOG_MAX_LEVEL) ;\
//     else fileLogMatrixPartial(level, mat, n, nS, nE, mS, mE, txt);
//
// #define FILE_LOGArray(level, arr, n, txt) \
//     if (level > FILELOG_MAX_LEVEL) ;\
//     else fileLogArray(level, arr, n, txt);
//
// #define FILE_LOGArrayPartial(level, arr, nS, nE, txt) \
//     if (level > FILELOG_MAX_LEVEL) ;\
//     else fileLogArrayPartial(level, arr, nS, nE, txt);
//
