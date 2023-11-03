#include "parallel_log.h"
#include "GlobalVariables.h"

#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

//#include <sys/time.h>
//#include <sys/time.h> Does not work on Windows
#include <ctime>
// This is for Windows build
#ifdef _WIN32
#include "winsock.h"
#include <Windows.h>
#endif
#include <stdint.h>


#include "mpi.h"

MPI_File logFilePtr;
// Char buffer for simple text log
char s1[256];
char s2[256];
va_list arg;

int gettimeofday(struct timeval* tp, struct timezone* tzp)
{
    // Note: some broken versions only have 8 trailing zero's, the correct epoch has 9 trailing zero's
    // This magic number is the number of 100 nanosecond intervals since January 1, 1601 (UTC)
    // until 00:00:00 January 1, 1970 
    static const uint64_t EPOCH = ((uint64_t)116444736000000000ULL);

/*    SYSTEMTIME  system_time;
    FILETIME    file_time;
    uint64_t    time;

    GetSystemTime(&system_time);
    SystemTimeToFileTime(&system_time, &file_time);
    time = ((uint64_t)file_time.dwLowDateTime);
    time += ((uint64_t)file_time.dwHighDateTime) << 32;

    tp->tv_sec = (long)((time - EPOCH) / 10000000L);
    tp->tv_usec = (long)(system_time.wMilliseconds * 1000);
    return 0;
    */
 #ifdef _WIN32 // Windows-specific code
    FILETIME file_time;
    GetSystemTimeAsFileTime(&file_time);
    
    ULARGE_INTEGER uli;
    uli.LowPart = file_time.dwLowDateTime;
    uli.HighPart = file_time.dwHighDateTime;
    uint64_t time = uli.QuadPart;
    
    struct timeval tv; // Renamed from tp to tv
    tv.tv_sec = (long)((time - EPOCH) / 10000000L);
    tv.tv_usec = (long)((time - EPOCH) % 10000000L);

    // Windows-specific code here

    #else // Code for non-Windows platforms (like Linux)
    struct timeval tv; // Renamed from tp to tv
    struct timespec ts;

    // Get the current time using clock_gettime on Linux
    clock_gettime(CLOCK_REALTIME, &ts);

    // Convert timespec to timeval
    tv.tv_sec = ts.tv_sec;
    tv.tv_usec = ts.tv_nsec / 1000;

    // Linux-specific code here

    #endif

    return 0;
}



void levelToString(enum logLevel level, char* result) {
    static const char* const buffer[] = { "ERROR  ", "WARNING", "INFO   ", "DEBUG  ", "DEBUG  " };
    char bufferTime[11];
    time_t t;
    time(&t);
    tm r;
    strftime(bufferTime, sizeof(buffer), "%X", localtime(&t));
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
