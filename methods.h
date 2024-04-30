#ifndef CONVERGENCE_METHODS_H
#define CONVERGENCE_METHODS_H

#include "auxfuncs.h"
#include "cuda_functions.h"

void runSimulation(char *method, real delta_t, real delta_x, real theta)
{
    // Number of steps
    int N = round(L / delta_x) + 1;               // Spatial steps (square tissue)
    int M = round(T / delta_t) + 1;                // Number of time steps

    // Allocate and populate time array
    real *time;
    time = (real *)malloc(M * sizeof(real));
    initializeTimeArray(time, M, delta_t);
    
    // Allocate and initialize the state variable
    real *V;
    V = (real *)malloc(N * N * sizeof(real));
    initializeStateVariable(V, N);

    // Auxiliary arrays for Thomas algorithm
    real *la = (real *)malloc(N * sizeof(real));    // subdiagonal
    real *lb = (real *)malloc(N * sizeof(real));    // diagonal
    real *lc = (real *)malloc(N * sizeof(real));    // superdiagonal

    // Populate auxiliary arrays for Thomas algorithm
    real D = sigma / (chi * Cm);
    real phi = (delta_t / (2 * delta_x * delta_x)) * D;       // For Thomas algorithm
    populateDiagonalThomasAlgorithm(la, lb, lc, N, phi);

    // Prefactorization
    thomasFactorConstantBatch(la, lb, lc, N);
    
    // Create directories
    char *pathToSaveData, *aux;
    pathToSaveData = (char *)malloc(MAX_STRING_SIZE * sizeof(char));
    aux = (char *)malloc(MAX_STRING_SIZE * sizeof(char));
    createDirectoriesAndFiles(method, theta, pathToSaveData, aux);

    // File names
    char infosFileName[MAX_STRING_SIZE];
    sprintf(infosFileName, "infos-%.8lf-%.6lf.txt", delta_t, delta_x);
    char lastFrameFileName[MAX_STRING_SIZE];
    sprintf(lastFrameFileName, "last-%.8lf-%.6lf.txt", delta_t, delta_x);

    // Infos file pointer
    FILE *fpInfos;
    sprintf(aux, "%s/%s", pathToSaveData, infosFileName);
    fpInfos = fopen(aux, "w"); 

    // CUDA variables and allocation
    real *d_V, *d_RHS, *d_Rv;
    real *d_la, *d_lb, *d_lc;

    CUDA_CALL(cudaMalloc(&d_V, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_RHS, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_Rv, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_la, N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_lb, N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_lc, N * sizeof(real)));

    // Copy memory of V from host to device of the matrices
    CUDA_CALL(cudaMemcpy(d_V, V, N * N * sizeof(real), cudaMemcpyHostToDevice));

    // Copy memory of diagonals from host to device
    CUDA_CALL(cudaMemcpy(d_la, la, N * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_lb, lb, N * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_lc, lc, N * sizeof(real), cudaMemcpyHostToDevice));

    // Block and grid size
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0); // Assuming device 0
    int max_block_dim_x = prop.maxThreadsDim[0];

    // For parallel Thomas
    printf("N = %d\n", N);
    int numBlocks = N / 100;
    if (numBlocks == 0)
        numBlocks = 1;
    int blockSize = round(N / numBlocks) + 1;        
    if (blockSize % 32 != 0)
        blockSize = 32 * ((blockSize / 32) + 1);
    if (blockSize > max_block_dim_x)
    {
        blockSize = max_block_dim_x;
        numBlocks = ceil(N / blockSize);
    }
    
    // For other kernels
    int GRID_SIZE = ceil((N*N*1.0) / (BLOCK_SIZE*1.0));
    if (GRID_SIZE == 0)
        GRID_SIZE = 1;

    printf("For 1st Part and Transpose -> Grid size %d, Block size %d\n", GRID_SIZE, BLOCK_SIZE);
    printf("Total for 1st Part and Transpose: %d\n", GRID_SIZE*BLOCK_SIZE);
    printf("For Thomas Algorithm -> Grid size: %d, Block size: %d\n", numBlocks, blockSize);
    printf("Total for Thomas Algorithm: %d\n", numBlocks*blockSize);
    printf("Spatial discretization N = %d!\n", N);
    printf("(2D) N * N = %d!\n", N*N);
    printf("Time discretization M = %d!\n", M);

    int timeStepCounter = 0;
    real timeStep = 0.0;
    if (strcmp(method, "SSI-ADI") == 0)
    {
        while (timeStepCounter < M)
        {
            // Get time step
            timeStep = time[timeStepCounter];

            // Solve the reaction and forcing term part
            parallelRHSForcing_SSI<<<GRID_SIZE, BLOCK_SIZE>>>(d_V, d_Rv, N, timeStep, delta_t, delta_x);
            cudaDeviceSynchronize();

            // Prepare right side of Thomas algorithm with explicit diffusion on j
            // Call the kernel
            prepareRHS<<<GRID_SIZE, BLOCK_SIZE>>>(d_V, d_RHS, d_Rv, N, phi); 
            cudaDeviceSynchronize();

            // Call the transpose kernel
            transposeDiagonalCol<<<GRID_SIZE, BLOCK_SIZE>>>(d_RHS, d_V, N);
            cudaDeviceSynchronize();
            
            // 1st: Implicit y-axis diffusion (lines)
            // Call the kernel
            cuThomasConstantBatch<<<numBlocks, blockSize>>>(d_la, d_lb, d_lc, d_V, N);
            cudaDeviceSynchronize();

            // Call the kernel
            prepareRHS<<<GRID_SIZE, BLOCK_SIZE>>>(d_V, d_RHS, d_Rv, N, phi); 
            cudaDeviceSynchronize();

            // Call the transpose kernel
            transposeDiagonalCol<<<GRID_SIZE, BLOCK_SIZE>>>(d_RHS, d_V, N);
            cudaDeviceSynchronize();

            // 2nd: Implicit x-axis diffusion (columns)                
            // Call the kernel
            cuThomasConstantBatch<<<numBlocks, blockSize>>>(d_la, d_lb, d_lc, d_V, N);
            cudaDeviceSynchronize();

            // Update time step counter
            timeStepCounter++;
        }
    }
    
    //Copy memory of d_V from device to host of the matrices (2D arrays)
    CUDA_CALL(cudaMemcpy(V, d_V, N * N * sizeof(real), cudaMemcpyDeviceToHost));

    // Save last frame
    FILE *fpLast;
    sprintf(aux, "%s/%s", pathToSaveData, lastFrameFileName);
    fpLast = fopen(aux, "w");
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
            fprintf(fpLast, "%e ", V[i * N + j]);
        fprintf(fpLast, "\n");
    }
    fclose(fpLast);

    // Write infos to file
    fprintf(fpInfos, "\nFor 1st Part and Transpose -> Grid size %d, Block size %d\n", GRID_SIZE, BLOCK_SIZE);
    fprintf(fpInfos, "Total threads: %d\n", GRID_SIZE*BLOCK_SIZE);
    ////////////////////////////////////////////////////
    fprintf(fpInfos, "\nFor 2nd Part -> Grid size: %d, Block size: %d\n", numBlocks, blockSize);
    fprintf(fpInfos, "Total threads: %d\n", numBlocks*blockSize);
    ////////////////////////////////////////////////////
    fprintf(fpInfos, "\ntheta = %lf\n", theta);
    fprintf(fpInfos, "L = %lf, T = %lf, N = %d, N*N = %d, M = %d\n", L, T, N, N*N, M);

    // Close files
    fclose(fpInfos);

    // Free memory
    free(time);

    // Free memory from host
    free(V);
    free(la);
    free(lb);
    free(lc);
    free(pathToSaveData);
    free(aux);

    // Free memory from device
    CUDA_CALL(cudaFree(d_V));
    CUDA_CALL(cudaFree(d_RHS));
    CUDA_CALL(cudaFree(d_Rv));
    CUDA_CALL(cudaFree(d_la));
    CUDA_CALL(cudaFree(d_lb));
    CUDA_CALL(cudaFree(d_lc));

    // Reset device
    CUDA_CALL(cudaDeviceReset());

    printf("Simulation finished!\n");

    return;
}


#endif // CONVERGENCE_METHODS_H