#ifndef INCLUDE_H
#define INCLUDE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <omp.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#define MAX_STRING_SIZE 100

#define BDIMX 16
#define BDIMY 16
#define BLOCK_SIZE 32

// Define real type via compile command line (-D{OPTION}, USE_DOUBLE or USE_FLOAT)
#ifdef USE_DOUBLE
typedef double real;
#define REAL_TYPE "double"
#else
typedef float real;
#define REAL_TYPE "float"
#endif

#define CUDA_CALL(call) \
do { \
    cudaError_t error = call; \
    if (error != cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d - %s\n", __FILE__, __LINE__, cudaGetErrorString(error)); \
        exit(EXIT_FAILURE); \
    } \
} while(0)


// Define problem via compile command line (-D{OPTION}):
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

// Define execution type via compile command line (-D{OPTION}, SERIAL or GPU)
#ifdef SERIAL
const real _pi = 3.14159265358979323846;
const int L = 1; // space

#ifdef LINMONO
const real T = 0.1;
// const real G = 1.5;         // omega^-1 * cm^-2
// const real sigma = 1.2e-3;  // omega^-1 * cm^-1
// const real chi = 1.0e3;     // cm^-1
// const real Cm = 1.0e-3;     // mF * cm^-2
const real G = 1.0;         // omega^-1 * cm^-2
const real sigma = 1.0;     // omega^-1 * cm^-1
const real chi = 1.0;       // cm^-1
const real Cm = 1.0;        // mF * cm^-2
#endif // LINMONO
#ifdef DIFFREAC
const real sigma = 1.0;
#endif // DIFFREAC
#ifdef DIFF
const real T = 0.5; // time
const real sigma = 1.0;
#endif // DIFF
#ifdef MONOAFHN
const real T = 0.1;
// const real G = 1.5;         // omega^-1 * cm^-2
// const real eta1 = 4.4;      // omega^-1 * cm^-1
// const real eta2 = 0.012;    // dimensionless
// const real eta3 = 1.0;      // dimensionless
// const real vth = 13.0;      // mV
// const real vp = 100.0;      // mV
// const real sigma = 1.2e-3;  // omega^-1 * cm^-1
// const real chi = 1.0e3;     // cm^-1
// const real Cm = 1.0e-3;     // mF * cm^-2
const real G = 1.0;         // omega^-1 * cm^-2
const real eta1 = 1.0;      // omega^-1 * cm^-1
const real eta2 = 1.0;      // dimensionless
const real eta3 = 1.0;      // dimensionless
const real vth = 1.0;       // mV
const real vp = 1.0;        // mV
const real sigma = 1.0;     // omega^-1 * cm^-1
const real chi = 1.0;       // cm^-1
const real Cm = 1.0;        // mF * cm^-2
#endif // MONOAFHN
#endif // SERIAL

#ifdef GPU
const __constant__ real _pi = 3.14159265358979323846;
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
const real T = 0.1;
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
#endif // GPU


#endif // INCLUDE_H