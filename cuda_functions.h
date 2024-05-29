#ifndef CUDA_FUNCTIONS_H
#define CUDA_FUNCTIONS_H

#include "include.h"

const __constant__ int d_L = 1;
const __constant__ real d_pi = 3.14159265358979323846;

#ifdef LINMONO
// __constant__ real d_G = 1.5;
// __constant__ real d_sigma = 1.2e-3; // omega^-1 * cm^-1
// __constant__ real d_chi = 1.0e3; // cm^-1
// __constant__ real d_Cm = 1.0e-3; // mF * cm^-2
__constant__ real d_G = 1.0;
__constant__ real d_sigma = 1.0; // omega^-1 * cm^-1
__constant__ real d_chi = 1.0; // cm^-1
__constant__ real d_Cm = 1.0; // mF * cm^-2
#endif // LINMONO

#ifdef DIFF
__constant__ real d_sigma = 1.0;
#endif // DIFF

#ifdef LINMONO
__host__ __device__ real exactSolution(real t, real x, real y)
{
    return (1.0-exp(-t)) * cos(d_pi*x/d_L) * cos(d_pi*y/d_L);
    // return (exp(-t)) * cos(d_pi*x) * cos(d_pi*y);
}

__host__ __device__ real forcingTerm(real x, real y, real t)
{
    // return cos(d_pi*x/d_L) * cos(d_pi*y/d_L) * (d_chi*d_Cm*exp(-t) + ((2.0*d_pi*d_pi*d_sigma)/(d_L*d_L))*(1.0-exp(-t)) + (d_chi*d_G)*(1.0-exp(-t)));
    // return exactSolution(t, x, y) * (d_chi*d_Cm + d_chi*d_G + 2.0*d_sigma*d_pi*d_pi/(d_L*d_L));
    return exactSolution(t, x, y) * (-1.0 + d_chi*d_G + 2.0*d_sigma*d_pi*d_pi/(d_L*d_L));
}
#endif // LINMONO

#ifdef DIFF
__host__ __device__ real exactSolution(real t, real x, real y)
{
    return (exp(-t)) * cos(d_pi*x) * cos(d_pi*y);
    // return exp(x+y-t); // Tese Ricardo
}
__host__ __device__ real forcingTerm(real x, real y, real t)
{
    // return cos(d_pi*x/d_L) * cos(d_pi*y/d_L) * (exp(-t) + ((2.0*d_pi*d_pi*d_sigma)/(d_L*d_L))*(1.0-exp(-t)));
    return exactSolution(t, x, y) * (-1.0+2.0*(d_pi*d_pi));
    // return exp(x+y-t) * (-1.0 - 2.0*d_sigma); // Tese Ricardo
}
#endif // DIFF

#ifdef PARALLEL
__global__ void initializeVariable(real *d_V, int N, real delta_x)
{
    unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;

    if (index < N*N)
    {
        unsigned int i = index / N;
        unsigned int j = index % N;

        real x = j * delta_x;
        real y = i * delta_x;

        d_V[index] = exactSolution(0.0, x, y);
    }
}

__global__ void errorXerror(real *d_V, real *d_out, int N, real t, real delta_x)
{
    unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;

    if (index < N*N)
    {
        unsigned int i = index / N;
        unsigned int j = index % N;

        real x = j * delta_x;
        real y = i * delta_x;

        real solution = exactSolution(t, x, y);

        d_out[index] = (d_V[index] - solution) * (d_V[index] - solution);
    }
}

__global__ void parallelRHSForcing_SSI(real *d_V, real *d_Rv, int N, real time, real delta_t, real delta_x)
{
    unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;

    if (index < N*N)
    {
        unsigned int i = index / N;
        unsigned int j = index % N;
        
        real actualV = d_V[index];

        // Diffusion component
        real diffusion = 0.0;

        // Diffusion in x
        if (j > 0 && j < N-1)
            diffusion = (d_V[i*N + j+1] - 2*actualV + d_V[i*N + j-1]) / (delta_x * delta_x);
        else if (j == 0)
            diffusion = (2*d_V[i*N + j+1] - 2*actualV) / (delta_x * delta_x);
        else if (j == N-1)
            diffusion = (2*d_V[i*N + j-1] - 2*actualV) / (delta_x * delta_x);

        // Diffusion in y
        if (i > 0 && i < N-1)
            diffusion += (d_V[(i+1)*N + j] - 2*actualV + d_V[(i-1)*N + j]) / (delta_x * delta_x);
        else if (i == 0)
            diffusion += (2*d_V[(i+1)*N + j] - 2*actualV) / (delta_x * delta_x);
        else if (i == N-1)
            diffusion += (2*d_V[(i-1)*N + j] - 2*actualV) / (delta_x * delta_x);

        // Forcing term at time t+(dt/2)
        real x = j * delta_x;
        real y = i * delta_x;
        real forcing = forcingTerm(x, y, time);

        // Aux variables
        real Vtilde = 0.0;
        real tildeRHS = 0.0;

        #ifdef LINMONO        
        // Calculate an approx for V
        Vtilde = actualV + (0.5 * delta_t) * (((d_sigma/(d_chi*d_Cm)) * diffusion) + (forcing/(d_chi*d_Cm)) - (d_G*actualV/d_Cm));

        // Now, recalculate the RHS with the approx
        tildeRHS = -(d_G*Vtilde/d_Cm);
        #endif // LINMONO

        // Update V reaction term
        d_Rv[index] = tildeRHS;
    }
}

__global__ void prepareRHS_diff_i(real *d_V, real *d_RHS, real *d_Rv, int N, real phi, real delta_t, real time, real delta_x)
{
    unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;

    if (index < N * N)
    {
        unsigned int i = index / N;
        unsigned int j = index % N;

        real actualV = d_V[index];
        real diffusion = 0.0;

        // Diffusion in y
        if (i > 0 && i < N-1)
            diffusion = (d_V[(i+1)*N + j] - 2*actualV + d_V[(i-1)*N + j]);
        else if (i == 0)
            diffusion = (2*d_V[(i+1)*N + j] - 2*actualV);
        else if (i == N-1)
            diffusion = (2*d_V[(i-1)*N + j] - 2*actualV);
        
        real x = j * delta_x;
        real y = i * delta_x;
        real forcing = forcingTerm(x, y, time); 

        #ifdef LINMONO
        d_RHS[index] = actualV + (phi * diffusion) + (0.5*delta_t*d_Rv[index]) + (0.5*delta_t*forcing/(d_chi*d_Cm));
        #endif // LINMONO
        #ifdef DIFF
        d_RHS[index] = actualV + (phi * diffusion) + (0.5*delta_t*forcing);
        #endif // DIFF
    }
}

__global__ void prepareRHS_diff_j(real *d_V, real *d_RHS, real *d_Rv, int N, real phi, real delta_t, real time, real delta_x)
{
    unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;

    if (index < N * N)
    {
        unsigned int i = index / N;
        unsigned int j = index % N;

        real actualV = d_V[index];
        real diffusion = 0.0;

        // Diffusion in y
        if (j > 0 && j < N-1)
            diffusion = (d_V[i*N + (j+1)] - 2*actualV + d_V[i*N + (j-1)]);
        else if (j == 0)
            diffusion = (2*d_V[i*N + (j+1)] - 2*actualV);
        else if (j == N-1)
            diffusion = (2*d_V[i*N + (j-1)] - 2*actualV);

        real x = j * delta_x;
        real y = i * delta_x;
        real forcing = forcingTerm(x, y, time); 

        #ifdef LINMONO
        d_RHS[index] = actualV + (phi * diffusion) + (0.5*delta_t*d_Rv[index]) + (0.5*delta_t*forcing/(d_chi*d_Cm));
        #endif // LINMONO
        #ifdef DIFF
        d_RHS[index] = actualV + (phi * diffusion) + (0.5*delta_t*forcing);
        #endif // DIFF
    }
}

__global__ void transposeDiagonalCol(real *in, real *out, unsigned int n)
{
    unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n * n)
    {
        unsigned int ix = i / n;
        unsigned int iy = i % n;
        out[iy * n + ix] = in[ix * n + iy];
    }
}

__global__ void cuThomasConstantBatch(real* la, real* lb, real* lc, real* d, unsigned int n)
{
	int rowCurrent;
	int rowPrevious;

	int rowAhead;

	// set the current row
	rowCurrent = threadIdx.x + blockDim.x*blockIdx.x;

	int i = 0;

	if ( rowCurrent < n ) 
	{
		//----- Forward Sweep
		d[rowCurrent] = d[rowCurrent] / lb[i];

		#pragma unroll
		for (i = 1; i < n; ++i) {
			rowPrevious = rowCurrent;
			rowCurrent += n;

			d[rowCurrent] = (d[rowCurrent] - la[i]*d[rowPrevious]) / (lb[i]);
		}

		//----- Back Sub
		d[rowCurrent] = d[rowCurrent];

		#pragma unroll
		for (i = n - 2; i >= 0; --i) {
			rowAhead    = rowCurrent;
			rowCurrent -= n;

			d[rowCurrent] = d[rowCurrent] - lc[i] * d[rowAhead];
		}
	}
}

__global__ void solveExplicitly(real *d_V, real *d_Vaux, int N, real t, real delta_t, real delta_x)
{
    unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;

    if (index < N*N)
    {
        unsigned int i = index / N;
        unsigned int j = index % N;
        
        real actualV = d_V[index];
        real diffusion = 0.0;

        // Diffusion component
        // Diffusion in x
        if (j > 0 && j < N-1)
            diffusion = (d_V[i*N + j+1] - 2*actualV + d_V[i*N + j-1]) / (delta_x * delta_x);
        else if (j == 0)
            diffusion = (2*d_V[i*N + j+1] - 2*actualV) / (delta_x * delta_x);
        else if (j == N-1)
            diffusion = (2*d_V[i*N + j-1] - 2*actualV) / (delta_x * delta_x);

        // Diffusion in y
        if (i > 0 && i < N-1)
            diffusion += (d_V[(i+1)*N + j] - 2*actualV + d_V[(i-1)*N + j]) / (delta_x * delta_x);
        else if (i == 0)
            diffusion += (2*d_V[(i+1)*N + j] - 2*actualV) / (delta_x * delta_x);
        else if (i == N-1)
            diffusion += (2*d_V[(i-1)*N + j] - 2*actualV) / (delta_x * delta_x);

        // Forcing term
        real x = j * delta_x;
        real y = i * delta_x;
        real forcing = forcingTerm(x, y, t);
        
        #ifdef LINMONO        
        d_Vaux[index] = actualV + delta_t * (((d_sigma/(d_chi*d_Cm)) * diffusion) + (forcing/(d_chi*d_Cm)) - (d_G*actualV/d_Cm));
        #endif // LINMONO
        #ifdef DIFF
        d_Vaux[index] = actualV + delta_t * ((d_sigma * diffusion) + forcing);
        #endif // DIFF
    }
}
#endif // PARALLEL

#endif // CUDA_FUNCTIONS_H