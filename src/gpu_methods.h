#ifndef GPU_METHODS_H
#define GPU_METHODS_H

#include "gpu_functions.h"
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
    #ifdef MONOAFHN
    real *W = (real *)malloc(N * N * sizeof(real));
    #endif // MONOAFHN
    initialize2DVariableWithExactSolution(V, N, delta_x); // TODO: just for now until testing is done
    #ifdef MONOAFHN
    initialize2DVariableWithValue(W, N, 0.0);
    #endif // MONOAFHN

    // Auxiliary arrays for Thomas algorithm
    real *la = (real *)malloc(N * sizeof(real));    // subdiagonal
    real *lb = (real *)malloc(N * sizeof(real));    // diagonal
    real *lc = (real *)malloc(N * sizeof(real));    // superdiagonal

    // Populate auxiliary arrays for Thomas algorithm
    real phi = (delta_t / (delta_x * delta_x));
    if (strcmp(method, "ADI") == 0 || strcmp(method, "SSI-ADI") == 0)
    {
        populateDiagonalThomasAlgorithm(la, lb, lc, N, 0.5*phi);
    }
    else if (strcmp(method, "theta-ADI") == 0)
    {
        populateDiagonalThomasAlgorithm(la, lb, lc, N, theta*phi);
    }

    // Prefactorization
    thomasFactorConstantBatch(la, lb, lc, N);

    // Create device variables
    real *d_V, *d_RHS, *d_Vtilde, *d_partRHS, *d_W;
    real *d_la, *d_lb, *d_lc;

    // Allocate memory on device
    CUDA_CALL(cudaMalloc(&d_V, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_RHS, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_Vtilde, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_partRHS, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_W, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_la, N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_lb, N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_lc, N * sizeof(real)));

    // Copy memory from host to device
    CUDA_CALL(cudaMemcpy(d_V, V, N * N * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_W, V, N * N * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_la, la, N * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_lb, lb, N * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_lc, lc, N * sizeof(real), cudaMemcpyHostToDevice));

    // CUDA variables and allocation
    int numBlocks = N / 100; 
    (numBlocks == 0) ? numBlocks = 1 : numBlocks;
    int blockSize = round(N / numBlocks) + 1;
    if (blockSize % 32 != 0)
        blockSize = 32 * ((blockSize / 32) + 1);
    
    // For ODEs
    int GRID_SIZE = ceil((N*N*1.0) / (BLOCK_SIZE*1.0));
    
    // Create directories
    char *pathToSaveData = (char *)malloc(MAX_STRING_SIZE * sizeof(char));
    char *aux = (char *)malloc(MAX_STRING_SIZE * sizeof(char));
    createDirectoriesAndFiles(method, theta, pathToSaveData, aux);

    // File names
    char infosFileName[MAX_STRING_SIZE];
    sprintf(infosFileName, "infos-%.8lf-%.6lf.txt", delta_t, delta_x);
    char lastFrameFileName[MAX_STRING_SIZE];
    sprintf(lastFrameFileName, "last-%.8lf-%.6lf.txt", delta_t, delta_x);
    char exactFileName[MAX_STRING_SIZE];
    sprintf(exactFileName, "exact-%.8lf-%.6lf.txt", delta_t, delta_x);
    char errorsFileName[MAX_STRING_SIZE];
    sprintf(errorsFileName, "errors-%.8lf-%.6lf.txt", delta_t, delta_x);

    // Infos file pointer
    sprintf(aux, "%s/%s", pathToSaveData, infosFileName);
    FILE *fpInfos = fopen(aux, "w");

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
            #ifdef SERIAL
            // ================================================!
            //  Calcula V em n + 1/2 -> Resultado vai para RHS !
            // ================================================!
            real x, y;
            for (int j = 0; j < N; j++)
            {
                for (int i = 0; i < N; i++)
                {   
                    x = j * delta_x;
                    y = i * delta_x;
                    #if defined LINMONO || defined DIFF
                    d[i] = 0.5*phi * V[i][lim(j-1,N)] + (1-2*0.5*phi) * V[i][j] + 0.5*phi * V[i][lim(j+1,N)] + 0.5 * delta_t * forcingTerm(x, y, actualTime);
                    #endif // LINMONO || DIFF
                }
                
                tridiag(la, lb, lc, c_prime, d_prime, N, d, result);
                for (int i = 0; i < N; i++)
                {
                    RHS[i][j] = result[i];
                }
            }

            // ================================================!
            //  Calcula V em n + 1 -> Resultado vai para V     !
            // ================================================!
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    x = j * delta_x;
                    y = i * delta_x;
                    #if defined LINMONO || defined DIFF
                    d[j] = 0.5*phi * RHS[lim(i-1,N)][j] + (1-2*0.5*phi) * RHS[i][j] + 0.5*phi * RHS[lim(i+1,N)][j] + 0.5 * delta_t * forcingTerm(x, y, actualTime+delta_t);
                    #endif // LINMONO || DIFF
                }
                
                tridiag(la, lb, lc, c_prime, d_prime, N, d, result);
                for (int j = 0; j < N; j++)
                {
                    V[i][j] = result[j];
                }
            }
            #endif // SERIAL

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
            computeApproxSSI<<<GRID_SIZE, BLOCK_SIZE>>>(N, delta_t, phi, delta_x, actualTime, d_V, d_Vtilde, d_partRHS, d_W);
            cudaDeviceSynchronize();

            // ================================================!
            //  Calculate V at n+1/2 -> Result goes to RHS     !
            //  diffusion implicit in x and explicit in y      !
            // ================================================!
            prepareRHSwithjDiffSSI<<<GRID_SIZE, BLOCK_SIZE>>>(N, phi, d_V, d_RHS, d_partRHS);
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
            prepareRHSwithiDiffSSI<<<GRID_SIZE, BLOCK_SIZE>>>(N, phi, d_Vtilde, d_V, d_partRHS);
            cudaDeviceSynchronize();

            parallelThomas<<<numBlocks, blockSize>>>(N, d_V, d_la, d_lb, d_lc);
            cudaDeviceSynchronize();

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
            #ifdef SERIAL
            // ================================================!
            //  Calcula Approx.                                !
            // ================================================!
            real x, y;
            for (int i = 0; i < N; i++)
            {   
                for (int j = 0; j < N; j++)
                {
                    x = j * delta_x;
                    y = i * delta_x;

                    #ifdef LINMONO
                    real diff_term = (sigma/(chi*Cm))*phi*(V[i][lim(j-1,N)] + V[lim(i-1,N)][j] - 4*V[i][j] + V[i][lim(j+1,N)] + V[lim(i+1,N)][j]);
                    real for_term = forcingTerm(x, y, actualTime+(0.5*delta_t))/(chi*Cm);
                    real reac_term = G*V[i][j]/Cm;
                    Vtilde[i][j] = V[i][j] + diff_term + (delta_t*(for_term - reac_term));

                    // Preparing part of the RHS of the following linear systems
                    real reac_tilde_term = G*Vtilde[i][j]/Cm;
                    partRHS[i][j] = delta_t*(for_term - ((1.0-theta)*reac_term) - (theta*reac_tilde_term));
                    #endif // LINMONO

                    #ifdef MONOAFHN
                    real diff_term = (sigma/(chi*Cm))*phi*(V[i][lim(j-1,N)] + V[lim(i-1,N)][j] - 4*V[i][j] + V[i][lim(j+1,N)] + V[lim(i+1,N)][j]);
                    real for_term = forcingTerm(x, y, actualTime+(0.5*delta_t), W[i][j])/(chi*Cm);
                    real RHS_V_term = RHS_V(V[i][j], W[i][j]);
                    Vtilde[i][j] = V[i][j] + diff_term + (delta_t*(for_term - RHS_V_term));

                    // Preparing part of the RHS of the following linear systems
                    real RHS_Vtilde_term = RHS_V(Vtilde[i][j], W[i][j]);
                    partRHS[i][j] = delta_t*(for_term - ((1.0-theta)*RHS_V_term) - (theta*RHS_Vtilde_term));

                    // Update Wn+1
                    real RHS_W_term = RHS_W(V[i][j], W[i][j]);
                    real Wtilde = W[i][j] + (delta_t*RHS_W_term);
                    W[i][j] = W[i][j] + delta_t*RHS_W(Vtilde[i][j], Wtilde);
                    #endif // MONOAFHN
                }
            }
            
            // ================================================!
            //  Calcula V em n + 1/2 -> Resultado vai para RHS !
            // ================================================!
            for (int j = 0; j < N; j++)
            {
                for (int i = 0; i < N; i++)
                {   
                    x = j * delta_x;
                    y = i * delta_x;

                    real diff_term = (sigma/(chi*Cm))*(1.0-theta)*phi*(V[i][lim(j-1,N)] - 2*V[i][j] + V[i][lim(j+1,N)]);
                    d[i] = V[i][j] + diff_term + 0.5*partRHS[i][j];
                }
                
                tridiag(la, lb, lc, c_prime, d_prime, N, d, result);
                for (int i = 0; i < N; i++)
                {
                    RHS[i][j] = result[i];
                }
            }

            // ================================================!
            //  Calcula V em n + 1 -> Resultado vai para V     !
            // ================================================!
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    x = j * delta_x;
                    y = i * delta_x;

                    real diff_term = (sigma/(chi*Cm))*(1.0-theta)*phi*(RHS[lim(i-1,N)][j] - 2*RHS[i][j] + RHS[lim(i+1,N)][j]);
                    d[j] = RHS[i][j] + diff_term + 0.5*partRHS[i][j];
                }
                
                tridiag(la, lb, lc, c_prime, d_prime, N, d, result);
                for (int j = 0; j < N; j++)
                {
                    V[i][j] = result[j];
                }
            }
            #endif // SERIAL

            // Update time step counter
            timeStepCounter++;
        }
    }
    #endif // LINMONO || MONOAFHN

    real finishTime = omp_get_wtime();
    real elapsedTime = finishTime - startTime;

    #ifdef GPU
    // // Get (v - solution)²
    // errorXerror<<<GRID_SIZE, BLOCK_SIZE>>>(d_V, d_RHS, N, actualTime, delta_x);
    
    // // Allocate 2D array for variable
    // real *temp = (real *)malloc(N * N * sizeof(real));

    // // Copy memory of d_V from device to host of the matrices (2D arrays)
    // CUDA_CALL(cudaMemcpy(temp, d_RHS, N * N * sizeof(real), cudaMemcpyDeviceToHost));

    // // Calculate the sum of errors²
    // real sum = 0.0;
    // for (int i = 0; i < N; i++)
    //     for (int j = 0; j < N; j++)
    //         sum += temp[i*N+j];
    
    // // Calculate norm-2 error
    // real norm2error = sqrt(delta_x*delta_x*sum);

    #endif // GPU
    // Copy memory of d_V from device to host V
    CUDA_CALL(cudaMemcpy(V, d_V, N * N * sizeof(real), cudaMemcpyDeviceToHost));

    real **exact = (real **)malloc(N * sizeof(real *));
    for (int i = 0; i < N; i++)
    {
        exact[i] = (real *)malloc(N * sizeof(real));
    }
    real norm2error = calculateNorm2Error(V, exact, N, T, delta_x);

    // Save last frame
    FILE *fpLast;
    sprintf(aux, "%s/%s", pathToSaveData, lastFrameFileName);
    fpLast = fopen(aux, "w");
    FILE *fpExact;
    sprintf(aux, "%s/%s", pathToSaveData, exactFileName);
    fpExact = fopen(aux, "w");
    FILE *fpErrors;
    sprintf(aux, "%s/%s", pathToSaveData, errorsFileName);
    fpErrors = fopen(aux, "w");
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            #ifdef GPU
            // fprintf(fpLast, "%e ", temp[i * N + j]);
            #endif // GPU
            // #ifdef SERIAL
            int index = i * N + j;
            fprintf(fpLast, "%e ", V[index]);
            fprintf(fpExact, "%e ", exact[i][j]);
            fprintf(fpErrors, "%e ", abs(V[index] - exact[i][j]));
            // #endif // SERIAL
        }
        fprintf(fpLast, "\n");
        fprintf(fpExact, "\n");
        fprintf(fpErrors, "\n");
    }
    fclose(fpLast);
    fclose(fpExact);
    fclose(fpErrors);

    // Write infos to file
    fprintf(fpInfos, "Domain Length = %d, Time = %f\n", L, T);
    fprintf(fpInfos, "delta_x = %lf, Space steps N = %d, N*N = %d\n", delta_x, N, N*N);
    fprintf(fpInfos, "delta_t = %lf, Time steps = %d\n", delta_t, M);
    fprintf(fpInfos, "Method %s\n", method);
    #ifdef GPU
    fprintf(fpInfos, "\nFor Approx, ODEs and Transpose -> Grid size %d, Block size %d\n", GRID_SIZE, BLOCK_SIZE);
    fprintf(fpInfos, "Total threads: %d\n", GRID_SIZE*BLOCK_SIZE);
    fprintf(fpInfos, "\nFor Thomas -> Grid size: %d, Block size: %d\n", numBlocks, blockSize);
    fprintf(fpInfos, "Total threads: %d\n", numBlocks*blockSize);
    #endif // GPU
    fprintf(fpInfos, "\nNorm-2 Error = %lf\n", norm2error);
    fprintf(fpInfos, "\nTotal execution time = %lf\n", elapsedTime);

    // Close files
    fclose(fpInfos);

    // Free memory
    free(time);
    for (int i = 0; i < N; i++)
    {
        free(exact[i]);
    }  
    free(exact);

    // Free memory from host
    #ifdef SERIAL
    for (int i = 0; i < N; i++)
    {
        free(V[i]);
        free(Vtilde[i]);
        free(RHS[i]);
        free(partRHS[i]);
        free(exact[i]);
        #ifdef MONOAFHN
        free(W[i]);
        #endif // MONOAFHN
    }
    free(V);
    free(Vtilde);
    free(RHS);
    free(partRHS);
    free(exact);
    free(c_prime);
    free(d_prime);
    free(d);
    free(result);
    #ifdef MONOAFHN
    free(W);
    #endif // MONOAFHN
    #endif // SERIAL
    free(V);
    free(Vtilde);
    free(RHS);
    free(partRHS);
    #ifdef MONOAFHN
    free(W);
    #endif // MONOAFHN
    free(la);
    free(lb);
    free(lc);
    free(pathToSaveData);
    free(aux);

    // Free memory from device
    CUDA_CALL(cudaFree(d_V));
    CUDA_CALL(cudaFree(d_Vtilde));
    CUDA_CALL(cudaFree(d_RHS));
    CUDA_CALL(cudaFree(d_partRHS));
    #ifdef MONOAFHN
    CUDA_CALL(cudaFree(d_W));
    #endif // MONOAFHN
    CUDA_CALL(cudaFree(d_la));
    CUDA_CALL(cudaFree(d_lb));
    CUDA_CALL(cudaFree(d_lc));

    // Reset device
    CUDA_CALL(cudaDeviceReset());

    printf("Simulation finished!\n");

    return;
}

#endif // GPU_METHODS_H