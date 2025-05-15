#ifndef CORE_DEFINITIONS_H
#define CORE_DEFINITIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <omp.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <time.h>
#include "logger.h"

#ifdef USE_CUDA

#include <cuda_runtime.h>

#endif // USE_CUDA

// Define path separator based on OS
#ifdef _WIN32
#include <direct.h>
#define mkdir(path, mode) _mkdir(path) // On Windows, _mkdir ignores mode
#define PATH_SEPARATOR '\\'
#else
#include <unistd.h>
#define PATH_SEPARATOR '/'
#endif

// ANSI color codes
#define RESET "\033[0m"
#define BLUE "\033[34m"
#define GREEN "\033[32m"
#define YELLOW "\033[33m"
#define RED "\033[31m"

// Macros for printing messages
#define INFOMSG(msg, ...) printf(BLUE "[i] " RESET msg, ##__VA_ARGS__)
#define SUCCESSMSG(msg, ...) printf(GREEN "[+] " RESET msg, ##__VA_ARGS__)
#define WARNINGMSG(msg, ...) printf(YELLOW "[!] " RESET msg, ##__VA_ARGS__)
#define ERRORMSG(msg, ...) printf(RED "\n[x] " RESET msg "\n", ##__VA_ARGS__)
#define DEBUGMSG(msg, ...) printf(YELLOW "[.] " RESET msg, ##__VA_ARGS__)
#define LINEMSG() printf(YELLOW "[.] " RESET "Line number %d in file %s\n", __LINE__, __FILE__)

// Define maximum string size
#define MAX_STRING_SIZE 200

// Define block size for GPU
#ifdef USE_CUDA

#define FULL_DOMAIN_BLOCK_SIZE_X 16
#define FULL_DOMAIN_BLOCK_SIZE_Y 16
#define THOMAS_KERNEL_BLOCK_SIZE 32

#endif // USE_CUDA

#define MAX_SYS_SIZE 10002
#define NUMTHREADS 6

// Convert CM to UM
#define CM_TO_UM(x) ((int)(x * 1.0e4))

// Define real type
typedef double real;
#define REAL_TYPE "double"
#define STRTOREAL strtod
#define FSCANF_REAL "%le"
#define PRINTF_REAL "%lf"
#ifdef USE_FLOAT
typedef float real;
#define REAL_TYPE "float"
#define STRTOREAL strtof
#define FSCANF_REAL "%e"
#define PRINTF_REAL "%f"
#endif

#define _PI 3.14159265358979323846f

// Define measurement structure
typedef struct
{
    real elapsedExecutionTime;
    real elapsedTime1stPart;
    real elapsedTime2ndPart;
    real elapsedTime1stLS;
    real elapsedTime2ndLS;
    real elapsedSaveFramesTime;
    real elapsedMeasureVelocityTime;
    real stimVelocity;
} Measurement;

// Define stimulus structure
typedef struct
{
    real amplitude;
    real begin_time;
    real duration;
    struct
    {
        real min;
        real max;
    } x_range, y_range, x_discretized, y_discretized;
} Stimulus;

// Define CUDA error checking
#define CUDA_CALL(call)                                                                                                        \
    do                                                                                                                         \
    {                                                                                                                          \
        cudaError_t error = call;                                                                                              \
        if (error != cudaSuccess)                                                                                              \
        {                                                                                                                      \
            fprintf(stderr, RED "\n[x] " RESET "CUDA error at %s:%d - %s\n\n", __FILE__, __LINE__, cudaGetErrorString(error)); \
            exit(EXIT_FAILURE);                                                                                                \
        }                                                                                                                      \
    } while (0)

#endif // CORE_DEFINITIONS_H