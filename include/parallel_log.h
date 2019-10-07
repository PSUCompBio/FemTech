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
#endif //__PARALLEL_LOG_H__
