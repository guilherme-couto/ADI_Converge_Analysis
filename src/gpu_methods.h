#ifndef GPU_METHODS_H
#define GPU_METHODS_H

#include "auxfuncs.h"

void runSimulationGPU(char *method, real delta_t, real delta_x, real theta)
{
    // Number of steps
    int N = round(L / delta_x) + 1;               // Spatial steps (square tissue)
    int M = round(T / delta_t);                   // Number of time steps

    // Allocate and populate time array
    real *time = (real *)malloc(M * sizeof(real));
    initializeTimeArray(time, M, delta_t);
    
    // Allocate arrays for variables
    real *V, *Vtilde, *RHS, *partRHS;
    V = (real *)malloc(N * N * sizeof(real));
    Vtilde = (real *)malloc(N * N * sizeof(real));
    RHS = (real *)malloc(N * N * sizeof(real));
    partRHS = (real *)malloc(N * N * sizeof(real));
    
    #ifdef CONVERGENCE_ANALYSIS
    initialize2DVariableWithExactSolution(V, N, delta_x);
    #else
    initialize2DVariableWithValue(V, N, V0);
    #endif // CONVERGENCE_ANALYSIS

    #ifdef MONOAFHN
    real *W = (real *)malloc(N * N * sizeof(real));
    initialize2DVariableWithValue(W, N, W0);
    #endif // MONOAFHN

    #ifdef INIT_WITH_SPIRAL
    initialize2DVariableFromFile(V, N, "./spiral_files/lastV_0.0005_0.0005.txt", delta_x);
    initialize2DVariableFromFile(W, N, "./spiral_files/lastW_0.0005_0.0005.txt", delta_x);
    #endif // INIT_WITH_SPIRAL

    // Auxiliary arrays for Thomas algorithm
    real *la = (real *)malloc(N * sizeof(real));    // subdiagonal
    real *lb = (real *)malloc(N * sizeof(real));    // diagonal
    real *lc = (real *)malloc(N * sizeof(real));    // superdiagonal

    // Populate auxiliary arrays for Thomas algorithm
    real phi = (delta_t / (delta_x * delta_x));
    if (strcmp(method, "ADI") == 0 || strcmp(method, "SSI-ADI") == 0)
    {
        populateDiagonalThomasAlgorithm(la, lb, lc, N, 0.5*phi*(sigma/(Cm*chi)));
    }
    else if (strcmp(method, "theta-ADI") == 0)
    {
        populateDiagonalThomasAlgorithm(la, lb, lc, N, theta*phi*(sigma/(Cm*chi)));
    }

    // Prefactorization
    thomasFactorConstantBatch(la, lb, lc, N);

    // Create device variables
    real *d_V, *d_RHS, *d_Vtilde, *d_partRHS;
    real *d_la, *d_lb, *d_lc;

    // Allocate memory on device
    CUDA_CALL(cudaMalloc(&d_V, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_RHS, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_Vtilde, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_partRHS, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_la, N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_lb, N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_lc, N * sizeof(real)));

    // Copy memory from host to device
    CUDA_CALL(cudaMemcpy(d_V, V, N * N * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_la, la, N * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_lb, lb, N * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_lc, lc, N * sizeof(real), cudaMemcpyHostToDevice));
    
    #ifdef MONOAFHN
    real *d_W;
    CUDA_CALL(cudaMalloc(&d_W, N * N * sizeof(real)));
    CUDA_CALL(cudaMemcpy(d_W, W, N * N * sizeof(real), cudaMemcpyHostToDevice));
    #endif // MONOAFHN

    #ifndef CONVERGENCE_ANALYSIS
    #ifdef MONOAFHN
    // Allocate array for the stimuli
    Stimulus *stimuli = (Stimulus *)malloc(numberOfStimuli * sizeof(Stimulus));
    populateStimuli(stimuli, delta_x);

    // Allocate and copy array for the stimuli on device
    Stimulus *d_stimuli;
    CUDA_CALL(cudaMalloc(&d_stimuli, numberOfStimuli * sizeof(Stimulus)));
    CUDA_CALL(cudaMemcpy(d_stimuli, stimuli, numberOfStimuli * sizeof(Stimulus), cudaMemcpyHostToDevice));
    #endif // MONOAFHN
    #endif // not def CONVERGENCE_ANALYSIS

    // CUDA grid and block allocation
    // For Thomas algorithm kernel
    int numBlocks = N / 100;
    if (numBlocks == 0)
        numBlocks = 1;

    int blockSize = round(N / numBlocks) + 1;
    if (blockSize % 32 != 0)
        blockSize = 32 * ((blockSize / 32) + 1);
    
    // For ODEs and Transpose kernels
    int GRID_SIZE = ceil((N*N*1.0) / (BLOCK_SIZE*1.0));

    // Create directories
    char *pathToSaveData = (char *)malloc(MAX_STRING_SIZE * sizeof(char));
    createDirectories(method, theta, pathToSaveData);
    
    #ifdef SAVE_FRAMES
    // Save frames
    char framesPath[MAX_STRING_SIZE];
    FILE *fpFrames;
    int frameSaveRate;
    snprintf(framesPath, MAX_STRING_SIZE*sizeof(char), "%s/frames_%.4f_%.4f.txt", pathToSaveData, delta_t, delta_x);
    fpFrames = fopen(framesPath, "w");
    frameSaveRate = ceil(M / numberOfFrames);
    
    #endif // SAVE_FRAMES
    
    int timeStepCounter = 0;
    real actualTime = 0.0;

    // Measure total execution time
    real startTime = omp_get_wtime();

    if (strcmp(method, "ADI") == 0)
    {
        while (timeStepCounter < M)
        {
            // Get time step
            actualTime = time[timeStepCounter];

            // Update time step counter
            timeStepCounter++;
        }
    }

    #if defined LINMONO || defined MONOAFHN
    else if (strcmp(method, "SSI-ADI") == 0)
    {
        while (timeStepCounter < M)
        {
            // Get time step
            actualTime = time[timeStepCounter];
            
            // ================================================!
            //  Calculate Approx. and ODEs                     !
            // ================================================!
            #ifdef CONVERGENCE_ANALYSIS
            computeApproxSSI<<<GRID_SIZE, BLOCK_SIZE>>>(N, delta_t, phi, delta_x, actualTime, d_V, d_Vtilde, d_partRHS, d_W);
            #else
            computeApproxSSI<<<GRID_SIZE, BLOCK_SIZE>>>(N, delta_t, phi, delta_x, actualTime, d_V, d_Vtilde, d_partRHS, d_W, d_stimuli);
            #endif // CONVERGENCE_ANALYSIS
            cudaDeviceSynchronize();

            // ================================================!
            //  Calculate V at n+1/2 -> Result goes to RHS     !
            //  diffusion implicit in x and explicit in y      !
            // ================================================!
            prepareRHSwithjDiff<<<GRID_SIZE, BLOCK_SIZE>>>(N, 0.5*phi, d_V, d_RHS, d_partRHS);
            cudaDeviceSynchronize();

            parallelThomas<<<numBlocks, blockSize>>>(N, d_RHS, d_la, d_lb, d_lc);
            cudaDeviceSynchronize();

            // ================================================!
            //  Transpose d_RHS to d_Vtilde                    !
            // ================================================!
            parallelTranspose<<<GRID_SIZE, BLOCK_SIZE>>>(N, d_RHS, d_Vtilde);
            cudaDeviceSynchronize();

            // ================================================!
            //  Calculate V at n+1 -> Result goes to V         !
            //  diffusion implicit in y and explicit in x      !
            // ================================================!
            prepareRHSwithiDiff<<<GRID_SIZE, BLOCK_SIZE>>>(N, 0.5*phi, d_Vtilde, d_V, d_partRHS);
            cudaDeviceSynchronize();

            parallelThomas<<<numBlocks, blockSize>>>(N, d_V, d_la, d_lb, d_lc);
            cudaDeviceSynchronize();

            #ifdef SAVE_FRAMES
            // If save frames is true and time step is multiple of frame save rate
            if (timeStepCounter % frameSaveRate == 0)
            {
                // Copy memory of d_V from device to host V
                CUDA_CALL(cudaMemcpy(V, d_V, N * N * sizeof(real), cudaMemcpyDeviceToHost));
                cudaDeviceSynchronize();

                // Save frame
                saveFrame(fpFrames, actualTime, V, N);
                printf("Frame at time %lf ms saved to %s\n", actualTime, framesPath);
            }
            #endif // SAVE_FRAMES

            // Update time step counter
            timeStepCounter++;
        }
    }
    
    else if (strcmp(method, "theta-ADI") == 0)
    {
        while (timeStepCounter < M)
        {
            // Get time step
            actualTime = time[timeStepCounter];
            
            // ================================================!
            //  Calculate Approx. and ODEs                     !
            // ================================================!
            #ifdef CONVERGENCE_ANALYSIS
            computeApproxthetaADI<<<GRID_SIZE, BLOCK_SIZE>>>(N, delta_t, phi, theta, delta_x, actualTime, d_V, d_Vtilde, d_partRHS, d_W);
            #else
            computeApproxthetaADI<<<GRID_SIZE, BLOCK_SIZE>>>(N, delta_t, phi, theta, delta_x, actualTime, d_V, d_Vtilde, d_partRHS, d_W, d_stimuli);
            #endif // CONVERGENCE_ANALYSIS
            cudaDeviceSynchronize();

            // ================================================!
            //  Calculate V at n+1/2 -> Result goes to RHS     !
            //  diffusion implicit in x and explicit in y      !
            // ================================================!
            prepareRHSwithjDiff<<<GRID_SIZE, BLOCK_SIZE>>>(N, (1.0-theta)*phi, d_V, d_RHS, d_partRHS);
            cudaDeviceSynchronize();

            parallelThomas<<<numBlocks, blockSize>>>(N, d_RHS, d_la, d_lb, d_lc);
            cudaDeviceSynchronize();

            // ================================================!
            //  Transpose d_RHS to d_Vtilde                    !
            // ================================================!
            parallelTranspose<<<GRID_SIZE, BLOCK_SIZE>>>(N, d_RHS, d_Vtilde);
            cudaDeviceSynchronize();

            // ================================================!
            //  Calculate V at n+1 -> Result goes to V         !
            //  diffusion implicit in y and explicit in x      !
            // ================================================!
            prepareRHSwithiDiff<<<GRID_SIZE, BLOCK_SIZE>>>(N, (1.0-theta)*phi, d_Vtilde, d_V, d_partRHS);
            cudaDeviceSynchronize();

            parallelThomas<<<numBlocks, blockSize>>>(N, d_V, d_la, d_lb, d_lc);
            cudaDeviceSynchronize();

            #ifdef SAVE_FRAMES
            // If save frames is true and time step is multiple of frame save rate
            if (timeStepCounter % frameSaveRate == 0)
            {
                // Copy memory of d_V from device to host V
                CUDA_CALL(cudaMemcpy(V, d_V, N * N * sizeof(real), cudaMemcpyDeviceToHost));
                cudaDeviceSynchronize();

                // Save frame
                saveFrame(fpFrames, actualTime, V, N);
                printf("Frame at time %lf ms saved to %s\n", actualTime, framesPath);
            }
            #endif // SAVE_FRAMES

            // Update time step counter
            timeStepCounter++;
        }
    }
    #endif // LINMONO || MONOAFHN

    real finishTime = omp_get_wtime();
    real elapsedTime = finishTime - startTime;

    // Copy memory of d_V from device to host V
    CUDA_CALL(cudaMemcpy(V, d_V, N * N * sizeof(real), cudaMemcpyDeviceToHost));

    // Calculate exact solution
    #ifdef CONVERGENCE_ANALYSIS
    real **exact = (real **)malloc(N * sizeof(real *));
    for (int i = 0; i < N; i++)
    {
        exact[i] = (real *)malloc(N * sizeof(real));
    }
    real norm2error = calculateNorm2Error(V, exact, N, T, delta_x);
    #endif // CONVERGENCE_ANALYSIS

    // Write infos to file
    char infosFilePath[MAX_STRING_SIZE];
    snprintf(infosFilePath, MAX_STRING_SIZE*sizeof(char), "%s/infos_%.4f_%.4f.txt", pathToSaveData, delta_t, delta_x);
    FILE *fpInfos = fopen(infosFilePath, "w");
    printf("Infos saved to %s\n", infosFilePath);
    fprintf(fpInfos, "Domain Length = %d, Time = %f\n", L, T);
    fprintf(fpInfos, "delta_x = %lf, Space steps N = %d, N*N = %d\n", delta_x, N, N*N);
    fprintf(fpInfos, "delta_t = %lf, Time steps = %d\n", delta_t, M);
    fprintf(fpInfos, "Method %s\n", method);
    fprintf(fpInfos, "\nFor Approx, ODEs and Transpose -> Grid size %d, Block size %d\n", GRID_SIZE, BLOCK_SIZE);
    fprintf(fpInfos, "Total threads: %d\n", GRID_SIZE*BLOCK_SIZE);
    fprintf(fpInfos, "\nFor Thomas -> Grid size: %d, Block size: %d\n", numBlocks, blockSize);
    fprintf(fpInfos, "Total threads: %d\n", numBlocks*blockSize);
    fprintf(fpInfos, "\nTotal execution time = %lf\n", elapsedTime);
    #ifdef CONVERGENCE_ANALYSIS
    fprintf(fpInfos, "\nNorm-2 Error = %lf\n", norm2error);
    #endif // CONVERGENCE_ANALYSIS
    fclose(fpInfos);

    // Save last frame
    char lastFrameFilePath[MAX_STRING_SIZE];
    snprintf(lastFrameFilePath, MAX_STRING_SIZE*sizeof(char), "%s/last_%.4f_%.4f.txt", pathToSaveData, delta_t, delta_x);
    FILE *fpLast = fopen(lastFrameFilePath, "w");
    printf("Last frame saved to %s\n", lastFrameFilePath);
    #ifdef CONVERGENCE_ANALYSIS
    char exactFilePath[MAX_STRING_SIZE];
    snprintf(exactFilePath, MAX_STRING_SIZE*sizeof(char), "%s/exact_%.4f_%.4f.txt", pathToSaveData, delta_t, delta_x);
    FILE *fpExact = fopen(exactFilePath, "w");
    printf("Exact solution saved to %s\n", exactFilePath);
    char errorsFilePath[MAX_STRING_SIZE];
    snprintf(errorsFilePath, MAX_STRING_SIZE*sizeof(char), "%s/errors_%.4f_%.4f.txt", pathToSaveData, delta_t, delta_x);
    FILE *fpErrors = fopen(errorsFilePath, "w");
    printf("Errors saved to %s\n", errorsFilePath);
    #endif // CONVERGENCE_ANALYSIS
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            int index = i * N + j;
            fprintf(fpLast, "%e ", V[index]);
            #ifdef CONVERGENCE_ANALYSIS
            fprintf(fpExact, "%e ", exact[i][j]);
            fprintf(fpErrors, "%e ", abs(V[index] - exact[i][j]));
            #endif // CONVERGENCE_ANALYSIS
        }
        fprintf(fpLast, "\n");
        #ifdef CONVERGENCE_ANALYSIS
        fprintf(fpExact, "\n");
        fprintf(fpErrors, "\n");
        #endif // CONVERGENCE_ANALYSIS
    }
    fclose(fpLast);
    #ifdef CONVERGENCE_ANALYSIS
    fclose(fpExact);
    fclose(fpErrors);
    #endif // CONVERGENCE_ANALYSIS

    // Free memory
    free(time);
    #ifdef CONVERGENCE_ANALYSIS
    for (int i = 0; i < N; i++)
    {
        free(exact[i]);
    }  
    free(exact);
    #endif // CONVERGENCE_ANALYSIS

    // Free memory from host
    free(V);
    free(Vtilde);
    free(RHS);
    free(partRHS);
    free(la);
    free(lb);
    free(lc);
    free(pathToSaveData);
    #ifdef MONOAFHN
    free(W);
    free(stimuli);
    #endif // MONOAFHN

    // Free memory from device
    CUDA_CALL(cudaFree(d_V));
    CUDA_CALL(cudaFree(d_Vtilde));
    CUDA_CALL(cudaFree(d_RHS));
    CUDA_CALL(cudaFree(d_partRHS));
    CUDA_CALL(cudaFree(d_la));
    CUDA_CALL(cudaFree(d_lb));
    CUDA_CALL(cudaFree(d_lc));
    #ifdef MONOAFHN
    CUDA_CALL(cudaFree(d_W));
    CUDA_CALL(cudaFree(d_stimuli));
    #endif // MONOAFHN

    // Reset device
    CUDA_CALL(cudaDeviceReset());

    return;
}

#endif // GPU_METHODS_H