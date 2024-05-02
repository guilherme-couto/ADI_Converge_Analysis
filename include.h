#ifndef INCLUDE_H
#define INCLUDE_H

#define MAX_STRING_SIZE 100

// Sizes of shared memory tile
#define BLOCK_SIZE 32

// Define real type
typedef double real;
#define REAL_TYPE "double"

#define SERIAL

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <omp.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

int L = 1;
int T = 1;    
real G = 1.5;         // omega^-1 * cm^-2
real sigma = 1.2e-3;  // omega^-1 * cm^-1
real chi = 1.0e3;     // cm^-1
real Cm = 1.0e-3;     // mF * cm^-2
real V_init = 0.0;

real _pi = 3.14159265358979323846;

#define CUDA_CALL(call) \
do { \
    cudaError_t error = call; \
    if (error != cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d - %s\n", __FILE__, __LINE__, cudaGetErrorString(error)); \
        exit(EXIT_FAILURE); \
    } \
} while(0)


#endif // INCLUDE_H