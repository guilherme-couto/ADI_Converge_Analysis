#ifndef CORE_DEFINITIONS_H
#define CORE_DEFINITIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <omp.h>
#ifdef GPU
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#endif // GPU
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <time.h>

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
#define FULL_DOMAIN_BLOCK_SIZE_X 16
#define FULL_DOMAIN_BLOCK_SIZE_Y 16
#define THOMAS_KERNEL_BLOCK_SIZE 32
#define MAX_SYS_SIZE 10002
#define NUMTHREADS 6

// Convert CM to UM
#define CM_TO_UM(x) ((int)(x * 1.0e4))

// Define real type
#ifdef USE_DOUBLE
typedef double real;
#define REAL_TYPE "double"
#define STRTOREAL strtod
#define FSCANF_REAL "%le"
#define PRINTF_REAL "%lf"
#else
typedef float real;
#define REAL_TYPE "float"
#define STRTOREAL strtof
#define FSCANF_REAL "%e"
#define PRINTF_REAL "%f"
#endif

// Define execution type
#ifdef SERIAL
#define EXECUTION_TYPE "SERIAL"
#define RUNSIMULATION runSimulationSerial
#endif
#ifdef OPENMP
#define EXECUTION_TYPE "OPENMP"
#define RUNSIMULATION runSimulationOpenMP
#endif
#ifdef GPU
#define EXECUTION_TYPE "GPU"
#define RUNSIMULATION runSimulationGPU
#endif

// Define file extension
#ifdef SAVE_AS_VTK
#define FILE_EXTENSION "vtk"
#define SAVEFRAME(filePath, Vm, Nx, Ny, dx, dy) saveVTK2D(filePath, Vm, Nx, Ny, dx, dy)
#else
#define SAVE_AS_TXT
#define FILE_EXTENSION "txt"
#define SAVEFRAME(file, Vm, Nx, Ny, dx, dy) saveFrame(file, Vm, Nx, Ny)
#endif

// Define CUDA error checking
#ifdef GPU
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
#endif

#endif // CORE_DEFINITIONS_H