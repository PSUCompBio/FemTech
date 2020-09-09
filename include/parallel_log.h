#ifndef __PARALLEL_LOG_H__
#define __PARALLEL_LOG_H__

#include <string.h>
#include <stdio.h>
#include <stdarg.h>

#include "mpi.h"

#ifndef FILELOG_MAX_LEVEL
#ifdef DEBUG
#define FILELOG_MAX_LEVEL DEBUGLOG
#else
#define FILELOG_MAX_LEVEL INFO
#endif
#endif

enum logLevel{ERROR, WARNING, INFO, DEBUGLOG, DEBUGLOGIGNORE};

const char* levelToString(enum logLevel level);
void initLog(const char *outFileName);
void finaliseLog();
void fileLog(enum logLevel level, const char* fmt, ...);
void fileLogMaster(enum logLevel level, const char* fmt, ...);
void fileLogSingle(enum logLevel level, const char* fmt, ...);
void fileLogMatrix(enum logLevel level, const double* mat, const int n, \
    const int m, const char* fmt, ...);
void fileLogMatrixRM(enum logLevel level, const double* mat, const int n, \
    const int m, const char* fmt, ...);
void fileLogMatrixSingle(enum logLevel level, const double* mat, const int n, \
    const int m, const char* fmt, ...);
void fileLogMatrixRMSingle(enum logLevel level, const double* mat, const int n, \
    const int m, const char* fmt, ...);
void fileLogArray(enum logLevel level, const double* arr, const int n, \
    const char* fmt, ...);
void fileLogArrayMaster(enum logLevel level, const double* arr, const int n, \
    const char* fmt, ...);
void fileLogArraySingle(enum logLevel level, const double* arr, const int n, \
    const char* fmt, ...);
void fileLogArrayInt(enum logLevel level, const int* arr, const int n, \
    const char* fmt, ...);
void fileLogArrayIntMaster(enum logLevel level, const int* arr, const int n, \
    const char* fmt, ...);
void fileLogArrayIntSingle(enum logLevel level, const int* arr, const int n, \
    const char* fmt, ...);

#ifndef FILELOG_MAX_LEVEL
#ifdef DEBUG
#define FILELOG_MAX_LEVEL DEBUGLOG
#else
#define FILELOG_MAX_LEVEL INFO
#endif
#endif

#define FILE_LOG(level, ...) \
    if (level > FILELOG_MAX_LEVEL) ;\
    else fileLog(level, __VA_ARGS__);

#define FILE_LOG_MASTER(level, ...) \
    if (level > FILELOG_MAX_LEVEL) ;\
    else fileLogMaster(level, __VA_ARGS__);

#define FILE_LOG_SINGLE(level, ...) \
    if (level > FILELOG_MAX_LEVEL) ;\
    else fileLogSingle(level, __VA_ARGS__);

#define FILE_LOGMatrix(level, mat, n, m, ...) \
    if (level > FILELOG_MAX_LEVEL) ;\
    else fileLogMatrix(level, mat, n, m, __VA_ARGS__);

#define FILE_LOGMatrixRM(level, mat, n, m, ...) \
    if (level > FILELOG_MAX_LEVEL) ;\
    else fileLogMatrixRM(level, mat, n, m, __VA_ARGS__);

#define FILE_LOGMatrix_SINGLE(level, mat, n, m, ...) \
    if (level > FILELOG_MAX_LEVEL) ;\
    else fileLogMatrixSingle(level, mat, n, m, __VA_ARGS__);

#define FILE_LOGMatrixRM_SINGLE(level, mat, n, m, ...) \
    if (level > FILELOG_MAX_LEVEL) ;\
    else fileLogMatrixRMSingle(level, mat, n, m, __VA_ARGS__);

#define FILE_LOGArray(level, arr, n, ...) \
    if (level > FILELOG_MAX_LEVEL) ;\
    else fileLogArray(level, arr, n, __VA_ARGS__);

#define FILE_LOGArrayMaster(level, arr, n, ...) \
    if (level > FILELOG_MAX_LEVEL) ;\
    else fileLogArrayMaster(level, arr, n, __VA_ARGS__);

#define FILE_LOGArraySingle(level, arr, n, ...) \
    if (level > FILELOG_MAX_LEVEL) ;\
    else fileLogArraySingle(level, arr, n, __VA_ARGS__);

#define FILE_LOGArrayInt(level, arr, n, ...) \
    if (level > FILELOG_MAX_LEVEL) ;\
    else fileLogArrayInt(level, arr, n, __VA_ARGS__);

#define FILE_LOGArrayIntMaster(level, arr, n, ...) \
    if (level > FILELOG_MAX_LEVEL) ;\
    else fileLogArrayIntMaster(level, arr, n, __VA_ARGS__);

#define FILE_LOGArrayIntSingle(level, arr, n, ...) \
    if (level > FILELOG_MAX_LEVEL) ;\
    else fileLogArrayIntSingle(level, arr, n, __VA_ARGS__);

#endif //__PARALLEL_LOG_H__
