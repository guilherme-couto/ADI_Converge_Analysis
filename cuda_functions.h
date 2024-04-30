#ifndef CUDA_FUNCTIONS_H
#define CUDA_FUNCTIONS_H

#include "include.h"

__constant__ real d_L = 1.0;
__constant__ real d_pi = 3.14159265358979323846;
__constant__ real d_G = 1.5;
__constant__ real d_sigma = 1.2e-3; // omega^-1 * cm^-1
__constant__ real d_chi = 1.0e3; // cm^-1
__constant__ real d_Cm = 1.0e-3; // mF * cm^-2

__global__ void parallelRHSForcing_SSI(real *d_V, real *d_Rv, int N, real t, real delta_t, real delta_x)
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
            diffusion = (d_V[i*N + j+1] - 2*d_V[index] + d_V[i*N + j-1]) / (delta_x * delta_x);
        else if (j == 0)
            diffusion = (2*d_V[i*N + j+1] - 2*d_V[index]) / (delta_x * delta_x);
        else if (j == N-1)
            diffusion = (2*d_V[i*N + j-1] - 2*d_V[index]) / (delta_x * delta_x);

        // Diffusion in y
        if (i > 0 && i < N-1)
            diffusion += (d_V[(i+1)*N + j] - 2*d_V[index] + d_V[(i-1)*N + j]) / (delta_x * delta_x);
        else if (i == 0)
            diffusion += (2*d_V[(i+1)*N + j] - 2*d_V[index]) / (delta_x * delta_x);
        else if (i == N-1)
            diffusion += (2*d_V[(i-1)*N + j] - 2*d_V[index]) / (delta_x * delta_x);

        // Forcing term
        real forcingTerm = 0.0;
        real x = j * delta_x;
        real y = i * delta_x;
        forcingTerm = cos(d_pi*x/d_L) * cos(d_pi*y/d_L) * (d_chi*d_Cm*exp(-t) + ((2.0*d_pi*d_pi*d_sigma)/(d_L*d_L))*(1.0-exp(-t)) + (d_chi*d_G)*(1.0-exp(-t))); 
        
        // Calculate the RHS with actual values
        real actualRHS = (forcingTerm/(d_chi*d_Cm)) - (d_G*actualV/d_Cm);
        
        // Calculate an approx for V
        real Vtilde;
        Vtilde = actualV + (0.5 * delta_t) * (((d_sigma/(d_chi*d_Cm)) * diffusion) + actualRHS);

        // Now, recalculate the forcing term at time t+(dt/2) and the RHS with the new values
        real new_t = t + (0.5 * delta_t);
        forcingTerm = cos(d_pi*x/d_L) * cos(d_pi*y/d_L) * (d_chi*d_Cm*exp(-new_t) + ((2.0*d_pi*d_pi*d_sigma)/(d_L*d_L))*(1.0-exp(-new_t)) + (d_chi*d_G)*(1.0-exp(-new_t)));
        real tildeRHS = (forcingTerm/(d_chi*d_Cm)) - (d_G*Vtilde/d_Cm);

        // Update V reaction term
        d_Rv[index] = tildeRHS;
    }
}

__global__ void prepareRHS(real *d_V, real *d_RHS, real *d_Rv, int N, real phi, real delta_t)
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
            diffusion += (d_V[(i+1)*N + j] - 2*actualV + d_V[(i-1)*N + j]);
        else if (i == 0)
            diffusion += (2*d_V[(i+1)*N + j] - 2*actualV);
        else if (i == N-1)
            diffusion += (2*d_V[(i-1)*N + j] - 2*actualV);

        d_RHS[index] = actualV + (phi * diffusion) + (0.5*delta_t*d_Rv[index]);
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
            diffusion = (d_V[i*N + j+1] - 2*d_V[index] + d_V[i*N + j-1]) / (delta_x * delta_x);
        else if (j == 0)
            diffusion = (2*d_V[i*N + j+1] - 2*d_V[index]) / (delta_x * delta_x);
        else if (j == N-1)
            diffusion = (2*d_V[i*N + j-1] - 2*d_V[index]) / (delta_x * delta_x);

        // Diffusion in y
        if (i > 0 && i < N-1)
            diffusion += (d_V[(i+1)*N + j] - 2*d_V[index] + d_V[(i-1)*N + j]) / (delta_x * delta_x);
        else if (i == 0)
            diffusion += (2*d_V[(i+1)*N + j] - 2*d_V[index]) / (delta_x * delta_x);
        else if (i == N-1)
            diffusion += (2*d_V[(i-1)*N + j] - 2*d_V[index]) / (delta_x * delta_x);

        // Forcing term
        real forcingTerm = 0.0;
        real x = j * delta_x;
        real y = i * delta_x;
        forcingTerm = cos(d_pi*x/d_L) * cos(d_pi*y/d_L) * (d_chi*d_Cm*exp(-t) + ((2.0*d_pi*d_pi*d_sigma)/(d_L*d_L))*(1.0-exp(-t)) + (d_chi*d_G)*(1.0-exp(-t)));
        
        real actualRHS = (forcingTerm/(d_chi*d_Cm)) - (d_G*actualV/d_Cm);
        
        d_Vaux[index] = actualV + delta_t * (((d_sigma/(d_chi*d_Cm)) * diffusion) + actualRHS);
    }
}

#endif // CUDA_FUNCTIONS_H