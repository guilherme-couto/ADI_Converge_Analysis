#ifndef GPU_FUNCTIONS_H
#define GPU_FUNCTIONS_H

#include "simulation_config.h"

#ifdef DIFF
#ifdef CONVERGENCE_ANALYSIS
__device__ real exactSolution(real t, real x, real y)
{
    return (exp(-t)) * cos(_pi*x) * cos(_pi*y);
}
__device__ real forcingTerm(real x, real y, real t)
{
    return exactSolution(t, x, y) * (-1.0+2.0*(_pi*_pi));
}
#endif // CONVERGENCE_ANALYSIS
#endif // DIFF

#ifdef LINMONO
#ifdef CONVERGENCE_ANALYSIS
__device__ real exactSolution(real t, real x, real y)
{
    return (exp(-t)) * cos(_pi*x/L) * cos(_pi*y/L);
}

__device__ real forcingTerm(real x, real y, real t)
{
    return exactSolution(t, x, y) * (-(chi*Cm) + chi*G + 2.0*(sigma/(chi*Cm))*_pi*_pi/(L*L));
}
#endif // CONVERGENCE_ANALYSIS
#endif // LINMONO

#ifdef MONODOMAIN
#if defined(CONVERGENCE_ANALYSIS) && defined(AFHN)
__host__ __device__ real exactSolution(real t, real x, real y)
{
    return (exp(-t)) * cos(_pi*x/L) * cos(_pi*y/L);
}

#ifdef AFHN
__device__ real forcingTerm(real x, real y, real t, real W)
{
    real exactV = exactSolution(t, x, y);
    real reaction = (G*exactV*(1.0-(exactV/vth)) * (1.0-(exactV/vp))) + (eta1*exactV*W);
    return (exactV * (-(chi*Cm) + 2.0*(sigma/(chi*Cm))*_pi*_pi/(L*L))) + (chi*reaction);
}
#endif // AFHN

// Kernel to compute the approximate solution of the reaction-diffusion system using the SSI-ADI
__global__ void computeApproxSSI(int N, real delta_t, real phi, real delta_x, real actualTime, real *d_V, real *d_Vtilde, real *d_partRHS, real *d_W)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    if (index < N * N)
    {
        int i = index / N;
        int j = index % N;

        real x = j * delta_x;
        real y = i * delta_x;

        int index_im1 = (i-1) * N + j;
        int index_ip1 = (i+1) * N + j;
        int index_jm1 = i * N + (j-1);
        int index_jp1 = i * N + (j+1);

        (i-1 == -1) ? (index_im1 = index_ip1) : index_im1;
        (i+1 == N) ? (index_ip1 = index_im1) : index_ip1;
        (j-1 == -1) ? (index_jm1 = index_jp1) : index_jm1;
        (j+1 == N) ? (index_jp1 = index_jm1) : index_jp1;

        real actualV = d_V[index];
        real im1V = d_V[index_im1];
        real ip1V = d_V[index_ip1];
        real jm1V = d_V[index_jm1];
        real jp1V = d_V[index_jp1];
        real diff_term = (sigma/(chi*Cm))*0.5*phi*(jm1V + im1V - 4*actualV + jp1V + ip1V);

        real actualW = d_W[index];

        real RHS_V_term = ((G*actualV*(1.0-(actualV/vth)) * (1.0-(actualV/vp))) + (eta1*actualV*actualW))/(Cm*chi);
        real for_term = forcingTerm(x, y, actualTime+(0.5*delta_t), actualW)/(chi*Cm);
        d_Vtilde[index] = actualV + diff_term + (0.5*delta_t*(for_term - RHS_V_term));

        real actualVtilde = d_Vtilde[index];

        // Preparing part of the RHS of the following linear systems
        real RHS_Vtilde_term = (G*actualVtilde*(1.0-(actualVtilde/vth)) * (1.0-(actualVtilde/vp))) + (eta1*actualVtilde*actualW);
        d_partRHS[index] = delta_t*(for_term - RHS_Vtilde_term);

        // Update Wn+1
        real RHS_W_term = eta2*((actualV/vp)-(eta3*actualW));
        real Wtilde = actualW + (0.5*delta_t*RHS_W_term);
        d_W[index] = actualW + delta_t*(eta2*((actualVtilde/vp)-(eta3*Wtilde)));
    }
}

// Kernel to compute the approximate solution of the reaction-diffusion system using the theta-ADI
__global__ void computeApproxthetaADI(int N, real delta_t, real phi, real theta, real delta_x, real actualTime, real *d_V, real *d_Vtilde, real *d_partRHS, real *d_W)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    if (index < N * N)
    {
        int i = index / N;
        int j = index % N;

        real x = j * delta_x;
        real y = i * delta_x;

        int index_im1 = (i-1) * N + j;
        int index_ip1 = (i+1) * N + j;
        int index_jm1 = i * N + (j-1);
        int index_jp1 = i * N + (j+1);

        real actualV = d_V[index];
        real im1V = d_V[index_im1];
        real ip1V = d_V[index_ip1];
        real jm1V = d_V[index_jm1];
        real jp1V = d_V[index_jp1];
        real diff_term = (sigma/(chi*Cm))*phi*(jm1V + im1V - 4*actualV + jp1V + ip1V);

        real actualV = d_V[index];
        real actualW = d_W[index];

        real for_term = forcingTerm(x, y, actualTime+(0.5*delta_t), actualW)/(chi*Cm);
        real RHS_V_term = ((G*actualV*(1.0-(actualV/vth)) * (1.0-(actualV/vp))) + (eta1*actualV*actualW))/(Cm*chi);
        d_Vtilde[index] = actualV + diff_term + (delta_t*(for_term - RHS_V_term));

        real actualVtilde = d_Vtilde[index];

        // Preparing part of the RHS of the following linear systems
        real RHS_Vtilde_term = (G*actualVtilde*(1.0-(actualVtilde/vth)) * (1.0-(actualVtilde/vp))) + (eta1*actualVtilde*actualW);
        d_partRHS[index] = delta_t*(for_term - ((1.0-theta)*RHS_V_term) - (theta*RHS_Vtilde_term));

        // Update Wn+1
        real RHS_W_term = eta2*((actualV/vp)-(eta3*actualW));
        real Wtilde = actualW + (delta_t*RHS_W_term);
        d_W[index] = actualW + delta_t*(eta2*((actualVtilde/vp)-(eta3*Wtilde)));
    }
}
#else
// Kernel to compute the approximate solution of the reaction-diffusion system using the SSI-ADI
__global__ void computeApproxSSI(int N, real delta_t, real phi, real delta_x, real actualTime, real *d_V, real *d_Vtilde, real *d_partRHS, stateVariables* d_sV, Stimulus *d_stimuli)
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

        (i-1 == -1) ? (index_im1 = index_ip1) : index_im1;
        (i+1 == N) ? (index_ip1 = index_im1) : index_ip1;
        (j-1 == -1) ? (index_jm1 = index_jp1) : index_jm1;
        (j+1 == N) ? (index_jp1 = index_jm1) : index_jp1;

        real actualV = d_V[index];
        real im1V = d_V[index_im1];
        real ip1V = d_V[index_ip1];
        real jm1V = d_V[index_jm1];
        real jp1V = d_V[index_jp1];
        real diff_term = (sigma/(chi*Cm))*0.5*phi*(jm1V + im1V - 4*actualV + jp1V + ip1V);

        // State variables
        #ifdef AFHN
        real actualW = d_sV[index].W;
        real RHS_V_term = ((G*actualV*(1.0-(actualV/vth)) * (1.0-(actualV/vp))) + (eta1*actualV*actualW))/(Cm*chi);
        #endif // AFHN

        real stim = 0.0;
        for (int si = 0; si < numberOfStimuli; si++)
        {
            if (actualTime >= d_stimuli[si].begin && actualTime <= d_stimuli[si].begin + d_stimuli[si].duration && j >= d_stimuli[si].xMinDisc && j <= d_stimuli[si].xMaxDisc && i >= d_stimuli[si].yMinDisc && i <= d_stimuli[si].yMaxDisc)
            {
                stim = d_stimuli[si].strength;
                break;
            }
        }

        real actualVtilde = actualV + diff_term + (0.5*delta_t*(stim - RHS_V_term));
        d_Vtilde[index] = actualVtilde;

        // Preparing part of the RHS of the following linear systems
        #ifdef AFHN
        real RHS_Vtilde_term = (G*actualVtilde*(1.0-(actualVtilde/vth)) * (1.0-(actualVtilde/vp))) + (eta1*actualVtilde*actualW);
        #endif // AFHN
        d_partRHS[index] = delta_t*(stim - RHS_Vtilde_term);

        // Update state variables
        #ifdef AFHN
        real RHS_W_term = eta2*((actualV/vp)-(eta3*actualW));
        real Wtilde = actualW + (0.5*delta_t*RHS_W_term);
        d_sV[index].W = actualW + delta_t*(eta2*((actualVtilde/vp)-(eta3*Wtilde)));
        #endif // AFHN
    }
}

// Kernel to compute the approximate solution of the reaction-diffusion system using the theta-ADI
__global__ void computeApproxthetaADI(int N, real delta_t, real phi, real theta, real delta_x, real actualTime, real *d_V, real *d_Vtilde, real *d_partRHS, stateVariables* d_sV, Stimulus *d_stimuli)
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

        (i-1 == -1) ? (index_im1 = index_ip1) : index_im1;
        (i+1 == N) ? (index_ip1 = index_im1) : index_ip1;
        (j-1 == -1) ? (index_jm1 = index_jp1) : index_jm1;
        (j+1 == N) ? (index_jp1 = index_jm1) : index_jp1;

        real actualV = d_V[index];
        real im1V = d_V[index_im1];
        real ip1V = d_V[index_ip1];
        real jm1V = d_V[index_jm1];
        real jp1V = d_V[index_jp1];
        real diff_term = (sigma/(chi*Cm))*phi*(jm1V + im1V - 4*actualV + jp1V + ip1V);

        // State variables
        #ifdef AFHN
        real actualW = d_sV[index].W;
        real RHS_V_term = ((G*actualV*(1.0-(actualV/vth)) * (1.0-(actualV/vp))) + (eta1*actualV*actualW))/(Cm*chi);
        #endif // AFHN
        
        real stim = 0.0;
        for (int si = 0; si < numberOfStimuli; si++)
        {
            if (actualTime >= d_stimuli[si].begin && actualTime <= d_stimuli[si].begin + d_stimuli[si].duration && j >= d_stimuli[si].xMinDisc && j <= d_stimuli[si].xMaxDisc && i >= d_stimuli[si].yMinDisc && i <= d_stimuli[si].yMaxDisc)
            {
                stim = d_stimuli[si].strength;
                break;
            }
        }

        real actualVtilde = actualV + diff_term + (delta_t*(stim - RHS_V_term));
        d_Vtilde[index] = actualVtilde;

        // Preparing part of the RHS of the following linear systems
        #ifdef AFHN
        real RHS_Vtilde_term = (G*actualVtilde*(1.0-(actualVtilde/vth)) * (1.0-(actualVtilde/vp))) + (eta1*actualVtilde*actualW);
        #endif // AFHN
        d_partRHS[index] = delta_t*(stim - ((1.0-theta)*RHS_V_term) - (theta*RHS_Vtilde_term));

        // Update state variables
        #ifdef AFHN
        real RHS_W_term = eta2*((actualV/vp)-(eta3*actualW));
        real Wtilde = actualW + (delta_t*RHS_W_term);
        d_sV[index].W = actualW + delta_t*(eta2*((actualVtilde/vp)-(eta3*Wtilde)));
        #endif // AFHN
    }
}
#endif

__global__ void prepareRHSwithiDiff(int N, real tau, real *d_V, real *d_RHS, real *d_partRHS)
{
    int index = blockDim.x * blockIdx.x + threadIdx.x;

    if (index < N * N)
    {
        int i = index / N;
        int j = index % N;

        int index_im1 = (i-1) * N + j;
        int index_ip1 = (i+1) * N + j;

        (i-1 == -1) ? (index_im1 = index_ip1) : index_im1;
        (i+1 == N) ? (index_ip1 = index_im1) : index_ip1;

        real actualV = d_V[index];

        real diff_term = (sigma/(chi*Cm))*tau*(d_V[index_im1] - 2*actualV + d_V[index_ip1]);
        d_RHS[index] = actualV + diff_term + 0.5*d_partRHS[index];
    }
}

__global__ void prepareRHSwithjDiff(int N, real tau, real *d_V, real *d_RHS, real *d_partRHS)
{
    int index = blockDim.x * blockIdx.x + threadIdx.x;

    if (index < N * N)
    {
        int i = index / N;
        int j = index % N;

        int transposedIndex = j * N + i;

        int index_jm1 = i * N + (j-1);
        int index_jp1 = i * N + (j+1);

        (j-1 == -1) ? (index_jm1 = index_jp1) : index_jm1;
        (j+1 == N) ? (index_jp1 = index_jm1) : index_jp1;

        real actualV = d_V[index];

        real diff_term = (sigma/(chi*Cm))*tau*(d_V[index_jm1] - 2*actualV + d_V[index_jp1]);
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
#endif // MONODOMAIN

#endif // GPU_FUNCTIONS_H