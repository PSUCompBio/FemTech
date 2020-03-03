#include "parallel_log.h"
#include "GlobalVariables.h"

#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <sys/time.h>

#include "mpi.h"

MPI_File logFilePtr;
// Char buffer for simple text log
char s1[256];
char s2[256];
va_list arg;

void levelToString(enum logLevel level, char* result) {
	static const char* const buffer[] = {"ERROR  ", "WARNING", "INFO   ", "DEBUG  ", "DEBUG  "};
  char bufferTime[11];
  time_t t;
  time(&t);
  tm r;
  strftime(bufferTime, sizeof(buffer), "%X", localtime_r(&t, &r));
  struct timeval tv;
  gettimeofday(&tv, 0);
  sprintf(result, "-%s.%03ld %s", bufferTime, (long)tv.tv_usec / 1000, buffer[level]); 
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
  char result[30] = {0};
  levelToString(level, result);
  sprintf(s1, "%s[%3d]: %s\n", result, world_rank, fmt);
  va_start(arg, fmt);
  vsprintf(s2, s1, arg);
  va_end(arg);
  MPI_File_write_ordered(logFilePtr, s2, strlen(s2), MPI_CHAR, MPI_STATUS_IGNORE);
}

void fileLogMaster(enum logLevel level, const char* fmt, ...) {
  if (world_rank == 0) {
    char result[30] = {0};
    levelToString(level, result);
    sprintf(s1, "%s[%3d]: %s\n", result, world_rank, fmt);
    va_start(arg, fmt);
    vsprintf(s2, s1, arg);
    va_end(arg);
    MPI_File_write_shared(logFilePtr, s2, strlen(s2), MPI_CHAR, MPI_STATUS_IGNORE);
  }
}

void fileLogSingle(enum logLevel level, const char* fmt, ...) {
  char result[30] = {0};
  levelToString(level, result);
  sprintf(s1, "%s[%3d]: %s\n", result, world_rank, fmt);
  va_start(arg, fmt);
  vsprintf(s2, s1, arg);
  va_end(arg);
  MPI_File_write_shared(logFilePtr, s2, strlen(s2), MPI_CHAR, MPI_STATUS_IGNORE);
}

void fileLogMatrix(enum logLevel level, const double* mat, const int n, \
    const int m, const char* fmt, ...) {
  char result[30] = {0};
  levelToString(level, result);
  sprintf(s1, "%s[%3d]: %s\n", result, world_rank, fmt);
  va_start(arg, fmt);
  vsprintf(s2, s1, arg);
  va_end(arg);
  int txtSize = strlen(s2)+(17*m+1)*n+10;
  char *sMat = (char*)malloc(txtSize*sizeof(char));
  int len = 0;
  len = sprintf(sMat, "%s", s2);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      len += sprintf(sMat+len, "%15.6E  ", mat[i+j*n]);
    }
    len += sprintf(sMat+len, "\n");
  }
  MPI_File_write_ordered(logFilePtr, sMat, strlen(sMat), MPI_CHAR, MPI_STATUS_IGNORE);
  free(sMat);
}

void fileLogMatrixRM(enum logLevel level, const double* mat, const int n, \
    const int m, const char* fmt, ...) {
  char result[30] = {0};
  levelToString(level, result);
  sprintf(s1, "%s[%3d]: %s\n", result, world_rank, fmt);
  va_start(arg, fmt);
  vsprintf(s2, s1, arg);
  va_end(arg);
  int txtSize = strlen(s2)+(17*m+1)*n+10;
  char *sMat = (char*)malloc(txtSize*sizeof(char));
  int len = 0;
  len = sprintf(sMat, "%s", s2);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      len += sprintf(sMat+len, "%15.6E  ", mat[j+i*m]);
    }
    len += sprintf(sMat+len, "\n");
  }
  MPI_File_write_ordered(logFilePtr, sMat, strlen(sMat), MPI_CHAR, MPI_STATUS_IGNORE);
  free(sMat);
}

void fileLogMatrixSingle(enum logLevel level, const double* mat, const int n, \
    const int m, const char* fmt, ...) {
  char result[30] = {0};
  levelToString(level, result);
  sprintf(s1, "%s[%3d]: %s\n", result, world_rank, fmt);
  va_start(arg, fmt);
  vsprintf(s2, s1, arg);
  va_end(arg);
  int txtSize = strlen(s2)+(17*m+1)*n+10;
  char *sMat = (char*)malloc(txtSize*sizeof(char));
  int len = 0;
  len = sprintf(sMat, "%s", s2);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      len += sprintf(sMat+len, "%15.6E  ", mat[i+j*n]);
    }
    len += sprintf(sMat+len, "\n");
  }
  MPI_File_write_shared(logFilePtr, sMat, strlen(sMat), MPI_CHAR, MPI_STATUS_IGNORE);
  free(sMat);
}

void fileLogMatrixRMSingle(enum logLevel level, const double* mat, const int n, \
    const int m, const char* fmt, ...) {
  char result[30] = {0};
  levelToString(level, result);
  sprintf(s1, "%s[%3d]: %s\n", result, world_rank, fmt);
  va_start(arg, fmt);
  vsprintf(s2, s1, arg);
  va_end(arg);
  int txtSize = strlen(s2)+(17*m+1)*n+10;
  char *sMat = (char*)malloc(txtSize*sizeof(char));
  int len = 0;
  len = sprintf(sMat, "%s", s2);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      len += sprintf(sMat+len, "%15.6E  ", mat[j+i*m]);
    }
    len += sprintf(sMat+len, "\n");
  }
  MPI_File_write_shared(logFilePtr, sMat, strlen(sMat), MPI_CHAR, MPI_STATUS_IGNORE);
  free(sMat);
}

void fileLogArray(enum logLevel level, const double* arr, const int n, \
    const char* fmt, ...) {
  char result[30] = {0};
  levelToString(level, result);
  sprintf(s1, "%s[%3d]: %s\n", result, world_rank, fmt);
  va_start(arg, fmt);
  vsprintf(s2, s1, arg);
  va_end(arg);
  int txtSize = strlen(s2)+(17*n+1)+10;
  char *sMat = (char*)malloc(txtSize*sizeof(char));
  int len = 0;
  len = sprintf(sMat, "%s", s2);
  for (int i = 0; i < n; ++i) {
    len += sprintf(sMat+len, "%15.6E  ", arr[i]);
  }
  len += sprintf(sMat+len, "\n");
  MPI_File_write_ordered(logFilePtr, sMat, strlen(sMat), MPI_CHAR, MPI_STATUS_IGNORE);
  free(sMat);
}

void fileLogArrayMaster(enum logLevel level, const double* arr, const int n, \
    const char* fmt, ...) {
  if (world_rank == 0) {
    char result[30] = {0};
    levelToString(level, result);
    sprintf(s1, "%s[%3d]: %s\n", result, world_rank, fmt);
    va_start(arg, fmt);
    vsprintf(s2, s1, arg);
    va_end(arg);
    int txtSize = strlen(s2)+(17*n+1)+10;
    char *sMat = (char*)malloc(txtSize*sizeof(char));
    int len = 0;
    len = sprintf(sMat, "%s", s2);
    for (int i = 0; i < n; ++i) {
      len += sprintf(sMat+len, "%15.6E  ", arr[i]);
    }
    len += sprintf(sMat+len, "\n");
    MPI_File_write_shared(logFilePtr, sMat, strlen(sMat), MPI_CHAR, MPI_STATUS_IGNORE);
    free(sMat);
  }
}

void fileLogArraySingle(enum logLevel level, const double* arr, const int n, \
    const char* fmt, ...) {
  char result[30] = {0};
  levelToString(level, result);
  sprintf(s1, "%s[%3d]: %s\n", result, world_rank, fmt);
  va_start(arg, fmt);
  vsprintf(s2, s1, arg);
  va_end(arg);
  int txtSize = strlen(s2)+(17*n+1)+10;
  char *sMat = (char*)malloc(txtSize*sizeof(char));
  int len = 0;
  len = sprintf(sMat, "%s", s2);
  for (int i = 0; i < n; ++i) {
    len += sprintf(sMat+len, "%15.6E  ", arr[i]);
  }
  len += sprintf(sMat+len, "\n");
  MPI_File_write_shared(logFilePtr, sMat, strlen(sMat), MPI_CHAR, MPI_STATUS_IGNORE);
  free(sMat);
}

void fileLogArrayInt(enum logLevel level, const int* arr, const int n, \
    const char* fmt, ...) {
  char result[30] = {0};
  levelToString(level, result);
  sprintf(s1, "%s[%3d]: %s\n", result, world_rank, fmt);
  va_start(arg, fmt);
  vsprintf(s2, s1, arg);
  va_end(arg);
  int txtSize = strlen(s2)+(13*n+1)+10;
  char *sMat = (char*)malloc(txtSize*sizeof(char));
  int len = 0;
  len = sprintf(sMat, "%s", s2);
  for (int i = 0; i < n; ++i) {
    len += sprintf(sMat+len, "%10d  ", arr[i]);
  }
  len += sprintf(sMat+len, "\n");
  MPI_File_write_ordered(logFilePtr, sMat, strlen(sMat), MPI_CHAR, MPI_STATUS_IGNORE);
  free(sMat);
}

void fileLogArrayIntMaster(enum logLevel level, const int* arr, const int n, \
    const char* fmt, ...) {
  if (world_rank == 0) {
    char result[30] = {0};
    levelToString(level, result);
    sprintf(s1, "%s[%3d]: %s\n", result, world_rank, fmt);
    va_start(arg, fmt);
    vsprintf(s2, s1, arg);
    va_end(arg);
    int txtSize = strlen(s2)+(13*n+1)+10;
    char *sMat = (char*)malloc(txtSize*sizeof(char));
    int len = 0;
    len = sprintf(sMat, "%s", s2);
    for (int i = 0; i < n; ++i) {
      len += sprintf(sMat+len, "%10d  ", arr[i]);
    }
    len += sprintf(sMat+len, "\n");
    MPI_File_write_shared(logFilePtr, sMat, strlen(sMat), MPI_CHAR, MPI_STATUS_IGNORE);
    free(sMat);
  }
}

void fileLogArrayIntSingle(enum logLevel level, const int* arr, const int n, \
    const char* fmt, ...) {
  char result[30] = {0};
  levelToString(level, result);
  sprintf(s1, "%s[%3d]: %s\n", result, world_rank, fmt);
  va_start(arg, fmt);
  vsprintf(s2, s1, arg);
  va_end(arg);
  int txtSize = strlen(s2)+(13*n+1)+10;
  char *sMat = (char*)malloc(txtSize*sizeof(char));
  int len = 0;
  len = sprintf(sMat, "%s", s2);
  for (int i = 0; i < n; ++i) {
    len += sprintf(sMat+len, "%10d  ", arr[i]);
  }
  len += sprintf(sMat+len, "\n");
  MPI_File_write_shared(logFilePtr, sMat, strlen(sMat), MPI_CHAR, MPI_STATUS_IGNORE);
  free(sMat);
}
