#ifndef GPU_FUNCTIONS_H
#define GPU_FUNCTIONS_H

#include "include.h"

#ifdef DIFF
__device__ real exactSolution(real t, real x, real y)
{
    return (exp(-t)) * cos(_pi*x) * cos(_pi*y);
}
__device__ real forcingTerm(real x, real y, real t)
{
    return exactSolution(t, x, y) * (-1.0+2.0*(_pi*_pi));
}
#endif // DIFF

#ifdef LINMONO
__device__ real exactSolution(real t, real x, real y)
{
    return (exp(-t)) * cos(_pi*x/L) * cos(_pi*y/L);
}

__device__ real forcingTerm(real x, real y, real t)
{
    return exactSolution(t, x, y) * (-(chi*Cm) + chi*G + 2.0*(sigma/(chi*Cm))*_pi*_pi/(L*L));
}
#endif // LINMONO

#ifdef MONOAFHN
__host__ __device__ real exactSolution(real t, real x, real y)
{
    return (exp(-t)) * cos(_pi*x/L) * cos(_pi*y/L);
}

__device__ real RHS_V(real V, real W)
{
    return (G*V*(1.0-(V/vth)) * (1.0-(V/vp))) + (eta1*V*W);
}

__device__ real RHS_W(real V, real W)
{
    return eta2*((V/vp)-(eta3*W));
}

__device__ real forcingTerm(real x, real y, real t, real W)
{
    real exactV = exactSolution(t, x, y);
    real reaction = (G*exactV*(1.0-(exactV/vth)) * (1.0-(exactV/vp))) + (eta1*exactV*W);
    return (exactV * (-(chi*Cm) + 2.0*(sigma/(chi*Cm))*_pi*_pi/(L*L))) + (chi*reaction);
}

// Kernel to compute the approximate solution of the reaction-diffusion system using the SSI-ADI
__global__ void computeApproxSSI(int N, real delta_t, real phi, real delta_x, real actualTime, real *d_V, real *d_Vtilde, real *d_partRHS, real *d_W)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    if (index < N * N)
    {
        int i = index / N;
        int j = index % N;

        int index_im1 = (i-1) * N + j;
        int index_ip1 = (i+1) * N + j;
        int index_jm1 = i * N + (j-1);
        int index_jp1 = i * N + (j+1);

        (i-1 == -1) ? index_im1 = index_ip1 : index_im1;
        (i+1 == N) ? index_ip1 = index_im1 : index_ip1;
        (j-1 == -1) ? index_jm1 = index_jp1 : index_jm1;
        (j+1 == N) ? index_jp1 = index_jm1 : index_jp1;

        real x = j * delta_x;
        real y = i * delta_x;

        real actualV = d_V[index];
        real actualW = d_W[index];

        real diff_term = (sigma/(chi*Cm))*0.5*phi*(d_V[index_jm1] + d_V[index_im1] - 4*actualV + d_V[index_jp1] + d_V[index_ip1]);
        real for_term = forcingTerm(x, y, actualTime+(0.5*delta_t), actualW)/(chi*Cm);
        real RHS_V_term = RHS_V(actualV, actualW);
        d_Vtilde[index] = actualV + diff_term + (0.5*delta_t*(for_term - RHS_V_term));

        real actualVtilde = d_Vtilde[index];

        // Preparing part of the RHS of the following linear systems
        real RHS_Vtilde_term = RHS_V(actualVtilde, actualW);
        d_partRHS[index] = delta_t*(for_term - RHS_Vtilde_term);

        // Update Wn+1
        real RHS_W_term = RHS_W(actualV, actualW);
        real Wtilde = actualW + (0.5*delta_t*RHS_W_term);
        d_W[index] = actualW + delta_t*RHS_W(actualVtilde, Wtilde);
    }
}

__global__ void prepareRHSwithiDiffSSI(int N, real phi, real *d_V, real *d_RHS, real *d_partRHS)
{
    int index = blockDim.x * blockIdx.x + threadIdx.x;

    if (index < N * N)
    {
        int i = index / N;
        int j = index % N;

        int index_im1 = (i-1) * N + j;
        int index_ip1 = (i+1) * N + j;

        (i-1 == -1) ? index_im1 = index_ip1 : index_im1;
        (i+1 == N) ? index_ip1 = index_im1 : index_ip1;

        real actualV = d_V[index];

        real diff_term = (sigma/(chi*Cm))*0.5*phi*(d_V[index_im1] - 2*actualV + d_V[index_ip1]);
        d_RHS[index] = actualV + diff_term + 0.5*d_partRHS[index];
    }
}

__global__ void prepareRHSwithjDiffSSI(int N, real phi, real *d_V, real *d_RHS, real *d_partRHS)
{
    int index = blockDim.x * blockIdx.x + threadIdx.x;

    if (index < N * N)
    {
        int i = index / N;
        int j = index % N;

        int transposedIndex = j * N + i;

        int index_jm1 = i * N + (j-1);
        int index_jp1 = i * N + (j+1);

        (j-1 == -1) ? index_jm1 = index_jp1 : index_jm1;
        (j+1 == N) ? index_jp1 = index_jm1 : index_jp1;

        real actualV = d_V[index];

        real diff_term = (sigma/(chi*Cm))*0.5*phi*(d_V[index_jm1] - 2*actualV + d_V[index_jp1]);
        d_RHS[transposedIndex] = actualV + diff_term + 0.5*d_partRHS[index];
    }
}

// From GLOSTER, Andrew et al. Efficient Interleaved Batch Matrix Solvers for CUDA. arXiv preprint arXiv:1909.04539, 2019.
// Kernel to solve the tridiagonal system using the Thomas algorithm
// N -> Number of rows
// d -> Result vector
// la -> Lower diagonal
// lb -> Main diagonal
// lc -> Upper diagonal
__global__ void parallelThomas(int N, real *d, real *la, real *lb, real *lc)
{
    int previousRow, nextRow;
    int currentRow = blockIdx.x * blockDim.x + threadIdx.x;
    int i = 0;

    if (currentRow < N)
    {
        // 1st: update auxiliary arrays
        d[currentRow] = d[currentRow] / lb[i];

#pragma unroll
        for (i = 1; i < N; i++)
        {
            previousRow = currentRow;
            currentRow += N;

            d[currentRow] = (d[currentRow] - la[i] * d[previousRow]) / (lb[i]);
        }

        // 2nd: update solution
        d[currentRow] = d[currentRow];

#pragma unroll
        for (i = N - 2; i >= 0; i--)
        {
            nextRow = currentRow;
            currentRow -= N;

            d[currentRow] = d[currentRow] - lc[i] * d[nextRow];
        }
    }
}

__global__ void parallelTranspose(int N, real *in, real *out)
{
    int index = blockDim.x * blockIdx.x + threadIdx.x;
    if (index < N * N)
    {
        int x = index / N;
        int y = index % N;
        out[y * N + x] = in[index];
    }
}
#endif // MONOAFHN

#endif // GPU_FUNCTIONS_H