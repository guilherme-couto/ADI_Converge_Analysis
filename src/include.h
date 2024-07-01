#ifndef INCLUDE_H
#define INCLUDE_H

#define MAX_STRING_SIZE 100

// Sizes of shared memory tile
#define BLOCK_SIZE 32

// Define real type
typedef double real;
#define REAL_TYPE "double"

// Define SERIAL or PARALLEL (PARALLEL is not working yet)
// #define PARALLEL // TODO: Implement parallel version
#define SERIAL

// Define problem:
// LINMONO -> Adapted monodomain with linear reaction (2D)
//            chi*Cm*dv/dt = sigma*Lap(v) - chi*G*v + forcing
//            Boundaries: Neumann
//
// DIFFREAC -> Diffusion with linear reaction (2D)
//             dv/dt + v = sigma*Lap(v) + forcing
//             Boundaries: Neumann
//
// DIFF -> Linear diffusion (2D)
//         dv/dt = sigma*Lap(v) + forcing
//         Boundaries: Neumann
//
// MONOAFHN -> Monodomain with adapted FitzHugh-Nagumo (2D)
//             { dv/dt = (sigma/(chi*Cm))*Lap(v) - RHS_v/Cm + forcing/(chi*Cm)
//             { dw/dt = RHS_w
//             RHS_v = (G*v*(1.0-(v/vth)) * (1.0-(v/vp))) + (eta1*v*w)
//             RHS_w = eta2 * ((v/vp)-(eta3*w))
//             Boundaries: Neumann
#define MONOAFHN

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <omp.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// Parameters
const __constant__ int L = 1; // space

#ifdef LINMONO
const __constant__ real T = 0.1;
// const __constant__ real G = 1.5;         // omega^-1 * cm^-2
// const __constant__ real sigma = 1.2e-3;  // omega^-1 * cm^-1
// const __constant__ real chi = 1.0e3;     // cm^-1
// const __constant__ real Cm = 1.0e-3;     // mF * cm^-2
const __constant__ real G = 1.0;         // omega^-1 * cm^-2
const __constant__ real sigma = 1.0;     // omega^-1 * cm^-1
const __constant__ real chi = 1.0;       // cm^-1
const __constant__ real Cm = 1.0;        // mF * cm^-2
#endif // LINMONO
#ifdef DIFFREAC
const __constant__ real sigma = 1.0;
#endif // DIFFREAC
#ifdef DIFF
const __constant__ real T = 0.5; // time
const __constant__ real sigma = 1.0;
#endif // DIFF
#ifdef MONOAFHN
const __constant__ real T = 0.1;
// const __constant__ real G = 1.5;         // omega^-1 * cm^-2
// const __constant__ real eta1 = 4.4;      // omega^-1 * cm^-1
// const __constant__ real eta2 = 0.012;    // dimensionless
// const __constant__ real eta3 = 1.0;      // dimensionless
// const __constant__ real vth = 13.0;      // mV
// const __constant__ real vp = 100.0;      // mV
// const __constant__ real sigma = 1.2e-3;  // omega^-1 * cm^-1
// const __constant__ real chi = 1.0e3;     // cm^-1
// const __constant__ real Cm = 1.0e-3;     // mF * cm^-2
const __constant__ real G = 1.0;         // omega^-1 * cm^-2
const __constant__ real eta1 = 1.0;      // omega^-1 * cm^-1
const __constant__ real eta2 = 1.0;      // dimensionless
const __constant__ real eta3 = 1.0;      // dimensionless
const __constant__ real vth = 1.0;       // mV
const __constant__ real vp = 1.0;        // mV
const __constant__ real sigma = 1.0;     // omega^-1 * cm^-1
const __constant__ real chi = 1.0;       // cm^-1
const __constant__ real Cm = 1.0;        // mF * cm^-2
#endif // MONOAFHN

const __constant__ real _pi = 3.14159265358979323846;

#define CUDA_CALL(call) \
do { \
    cudaError_t error = call; \
    if (error != cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d - %s\n", __FILE__, __LINE__, cudaGetErrorString(error)); \
        exit(EXIT_FAILURE); \
    } \
} while(0)


#endif // INCLUDE_H