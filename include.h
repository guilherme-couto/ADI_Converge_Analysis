#ifndef INCLUDE_H
#define INCLUDE_H

#define MAX_STRING_SIZE 100

// Sizes of shared memory tile
#define BLOCK_SIZE 32

// Define real type
typedef double real;
#define REAL_TYPE "double"

// Define SERIAL or PARALLEL
#define PARALLEL

// Define problem:
// LINMONO -> Adapted monodomain with linear reaction (2D)
//            chi*(Cm*dv/dt + G*v) = sigma*Lap(v) + forcing
//            v(x,0) = 0; dv(0,t)/dx = dv(L,t)/dx = 0 (Neumann)
//
// DIFFREAC -> Diffusion with linear reaction (2D)
//             dv/dt + v = sigma*Lap(v) + forcing
//             v(x,0) = 0; dv(0,t)/dx = dv(L,t)/dx = 0 (Neumann)
//
// DIFF -> Linear diffusion (2D)
//         dv/dt = sigma*Lap(v) + forcing
//         dv(0,t)/dx = dv(L,t)/dx = 0 (Neumann)
//
// SYSDIFFREAC -> System of equations with diffusion with reaction (2D)
//                du/ht - nu Lap(u) + (u^2 v) = forcing1(x,y,t)
//                dv/ht - nu Lap(v) + (u v^2) = forcing2(x,y,t)
//                v(x,0) = 0; dv(0,t)/dx = dv(L,t)/dx = 0 (Neumann)
//

#define DIFF

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <omp.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// Domain
int L = 1; // space
real T = 1; // time

// Parameters
#ifdef LINMONO
real G = 1.5;         // omega^-1 * cm^-2
real sigma = 1.2e-3;  // omega^-1 * cm^-1
real chi = 1.0e3;     // cm^-1
real Cm = 1.0e-3;     // mF * cm^-2
#endif // LINMONO
#ifdef DIFFREAC
real sigma = 1.0;
real V_init = 0.0;
#endif // DIFFREAC
#ifdef DIFF
real sigma = 1.0;
#endif // DIFF

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