#ifndef CONVERGENCE_METHODS_H
#define CONVERGENCE_METHODS_H

#include "cuda_functions.h"
#include "auxfuncs.h"

void runSimulation(char *method, real delta_t, real delta_x, real theta)
{
    // Number of steps
    int N = round(L / delta_x) + 1;               // Spatial steps (square tissue)
    int M = round(T / delta_t) + 1;               // Number of time steps

    // Allocate and populate time array
    real *time = (real *)malloc(M * sizeof(real));
    initializeTimeArray(time, M, delta_t);
    
    #ifdef SERIAL
    // Allocate 2D arrays for variables
    real **V, **Rv, **RHS, **exact;
    V = (real **)malloc(N * sizeof(real *));
    Rv = (real **)malloc(N * sizeof(real *));
    RHS = (real **)malloc(N * sizeof(real *));
    exact = (real **)malloc(N * sizeof(real *));
    real *c_prime = (real *)malloc(N * sizeof(real));   // aux for Thomas
    real *d_prime = (real *)malloc(N * sizeof(real));   // aux for Thomas
    real *d = (real *)malloc(N * sizeof(real));
    real *result = (real *)malloc(N * sizeof(real));
    for (int i = 0; i < N; i++)
    {
        V[i] = (real *)malloc(N * sizeof(real));
        Rv[i] = (real *)malloc(N * sizeof(real));
        RHS[i] = (real *)malloc(N * sizeof(real));
        exact[i] = (real *)malloc(N * sizeof(real));
    }
    initializeStateVariable(V, N, delta_x);
    #endif // SERIAL

    // Auxiliary arrays for Thomas algorithm
    real *la = (real *)malloc(N * sizeof(real));    // subdiagonal
    real *lb = (real *)malloc(N * sizeof(real));    // diagonal
    real *lc = (real *)malloc(N * sizeof(real));    // superdiagonal

    // Populate auxiliary arrays for Thomas algorithm
    #ifdef LINMONO
    real D = sigma / (chi * Cm);
    #endif // LINMONO
    #ifdef DIFF
    real D = sigma;
    #endif // DIFF
    real phi = (delta_t / (2 * delta_x * delta_x));       // For Thomas algorithm
    populateDiagonalThomasAlgorithm(la, lb, lc, N, phi);

    // Prefactorization
    #ifdef PARALLEL
    thomasFactorConstantBatch(la, lb, lc, N);
    #endif // PARALLELL
    
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

    // CUDA variables and allocation
    #ifdef PARALLEL
    real *d_V, *d_RHS, *d_Rv;
    real *d_la, *d_lb, *d_lc;

    CUDA_CALL(cudaMalloc(&d_V, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_RHS, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_Rv, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_la, N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_lb, N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_lc, N * sizeof(real)));

    // Copy memory of diagonals from host to device
    CUDA_CALL(cudaMemcpy(d_la, la, N * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_lb, lb, N * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_lc, lc, N * sizeof(real), cudaMemcpyHostToDevice));

    // Block and grid size
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0); // Assuming device 0
    int max_block_dim_x = prop.maxThreadsDim[0];

    // Blocks and threads for parallel Thomas (N calls)
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
    
    // Blocks and threads for kernels that go through all elements (N*N calls)
    int GRID_SIZE = ceil((N*N*1.0) / (BLOCK_SIZE*1.0));
    if (GRID_SIZE == 0)
        GRID_SIZE = 1;

    // Initialize state variable with exact solution at t = 0
    initializeVariable<<<GRID_SIZE, BLOCK_SIZE>>>(d_V, N, delta_x);
    #endif // PARALLEL

    int timeStepCounter = 0;
    real actualTime = 0.0;
    if (strcmp(method, "SSI-ADI") == 0)
    {
        while (timeStepCounter < M)
        {
            // Get time step
            actualTime = time[timeStepCounter];

            #ifdef PARALLEL
            // Solve the reaction and forcing term part
            parallelRHSForcing_SSI<<<GRID_SIZE, BLOCK_SIZE>>>(d_V, d_Rv, N, actualTime+(0.5*delta_t), delta_t, delta_x);
            cudaDeviceSynchronize();

            // Prepare right side of Thomas algorithm with explicit diffusion on y
            // Call the kernel
            prepareRHS_diff_i<<<GRID_SIZE, BLOCK_SIZE>>>(d_V, d_RHS, d_Rv, N, phi, delta_t, actualTime+(0.5*delta_t), delta_x); 
            cudaDeviceSynchronize();

            // Call the transpose kernel
            transposeDiagonalCol<<<GRID_SIZE, BLOCK_SIZE>>>(d_RHS, d_V, N);
            cudaDeviceSynchronize();
            
            // 1st: Implicit x-axis diffusion (lines)
            // Call the kernel
            cuThomasConstantBatch<<<numBlocks, blockSize>>>(d_la, d_lb, d_lc, d_V, N);
            cudaDeviceSynchronize();

            // Prepare right side of Thomas algorithm with explicit diffusion on x
            // Call the kernel
            prepareRHS_diff_j<<<GRID_SIZE, BLOCK_SIZE>>>(d_V, d_RHS, d_Rv, N, phi, delta_t, actualTime+(0.5*delta_t), delta_x); 
            cudaDeviceSynchronize();

            // 2nd: Implicit y-axis diffusion (columns)                
            // Call the kernel
            cuThomasConstantBatch<<<numBlocks, blockSize>>>(d_la, d_lb, d_lc, d_RHS, N);
            cudaDeviceSynchronize();

            // Copy d_RHS to d_V
            CUDA_CALL(cudaMemcpy(d_V, d_RHS, N * N * sizeof(real), cudaMemcpyDeviceToDevice));
            #endif // PARALLEL

            #ifdef SERIAL
            // Calculate the approx of V, the forcing term and reaction
            calculateVApprox(V, Rv, N, delta_x, delta_t, actualTime);

            // Prepare RHS for 1st part of ADI
            // RHS will have the contribution of the diffusion through y (columns of the matrix)
            prepareRHS_explicit_y(V, RHS, Rv, N, phi, delta_t, actualTime+(0.5*delta_t), delta_x);

            // Call Thomas
            // The solution of the system will give the diffusion through x (lines of the matrix)
            for (int i = 0; i < N; i++)
                thomasAlgorithm(la, lb, lc, c_prime, d_prime, N, RHS[i]);

            copyMatrices(RHS, V, N);
            
            // Prepare RHS for 2nd part of ADI
            // RHS will have the contribution of the diffusion through x (lines of the matrix)
            prepareRHS_explicit_x(V, RHS, Rv, N, phi, delta_t, actualTime+(0.5*delta_t), delta_x);

            // Call Thomas
            // The solution of the system will give the diffusion through y (columns of the matrix)
            for (int j = 0; j < N; j++)
            {
                copyColumnToVector(RHS, d, N, j);
                thomasAlgorithm(la, lb, lc, c_prime, d_prime, N, d);
                copyVectorToColumn(V, d, N, j);
            }
            #endif // SERIAL

            // Update time step counter
            timeStepCounter++;
        }
    }

    else if (strcmp(method, "FE") == 0)
    {
        while (timeStepCounter < M)
        {
            // Get time step
            actualTime = time[timeStepCounter];

            #ifdef PARALLEL
            // Solve the reaction and forcing term part
            solveExplicitly<<<GRID_SIZE, BLOCK_SIZE>>>(d_V, d_Rv, N, actualTime, delta_t, delta_x);
            cudaDeviceSynchronize();

            // Copy d_Rv to d_V
            CUDA_CALL(cudaMemcpy(d_V, d_Rv, N * N * sizeof(real), cudaMemcpyDeviceToDevice));
            #endif // PARALLEL
            
            #ifdef SERIAL
            if (actualTime < delta_t)
            {
                printf("%f/%e\n", actualTime, V[0][0]);
                printf("dx*dx = %e\n", delta_x*delta_x);
            }
            solveExplicitly(V, RHS, N, delta_x, delta_t, actualTime);
            #endif // SERIAL

            // Update time step counter
            timeStepCounter++;
        }
    }
    
    else if (strcmp(method, "ADI") == 0)
    {
        while (timeStepCounter < M)
        {
            // Get time step
            actualTime = time[timeStepCounter];

            #ifdef PARALLEL
            // Prepare right side of Thomas algorithm with explicit diffusion on y
            // Call the kernel
            prepareRHS_diff_i<<<GRID_SIZE, BLOCK_SIZE>>>(d_V, d_RHS, d_Rv, N, phi, delta_t, actualTime+(0.5*delta_t), delta_x); 
            cudaDeviceSynchronize();

            // Call the transpose kernel
            transposeDiagonalCol<<<GRID_SIZE, BLOCK_SIZE>>>(d_RHS, d_V, N);
            cudaDeviceSynchronize();
            
            // 1st: Implicit x-axis diffusion (lines)
            // Call the kernel
            cuThomasConstantBatch<<<numBlocks, blockSize>>>(d_la, d_lb, d_lc, d_V, N);
            cudaDeviceSynchronize();

            // Prepare right side of Thomas algorithm with explicit diffusion on x
            // Call the kernel
            prepareRHS_diff_j<<<GRID_SIZE, BLOCK_SIZE>>>(d_V, d_RHS, d_Rv, N, phi, delta_t, actualTime+(0.5*delta_t), delta_x); 
            cudaDeviceSynchronize();

            // 2nd: Implicit y-axis diffusion (columns)                
            // Call the kernel
            cuThomasConstantBatch<<<numBlocks, blockSize>>>(d_la, d_lb, d_lc, d_RHS, N);
            cudaDeviceSynchronize();

            // Copy d_RHS to d_V
            CUDA_CALL(cudaMemcpy(d_V, d_RHS, N * N * sizeof(real), cudaMemcpyDeviceToDevice));
            #endif // PARALLEL

            #ifdef SERIAL
            // // Prepare RHS for 1st part of ADI
            // // RHS will have the contribution of the diffusion through x (lines of the matrix)
            // prepareRHS_explicit_x(V, RHS, Rv, N, phi, delta_t, actualTime+(0.0*delta_t), delta_x);

            // // Call Thomas
            // // The solution of the system will give the diffusion through y (columns of the matrix)
            // for (int j = 0; j < N; j++)
            // {
            //     copyColumnToVector(RHS, d, N, j);
            //     thomasAlgorithm(la, lb, lc, c_prime, d_prime, N, d);
            //     copyVectorToColumn(V, d, N, j);
            // }

            // // Prepare RHS for 2nd part of ADI
            // // RHS will have the contribution of the diffusion through y (columns of the matrix)
            // prepareRHS_explicit_y(V, RHS, Rv, N, phi, delta_t, actualTime+(1.0*delta_t), delta_x);
            
            // // Call Thomas
            // // The solution of the system will give the diffusion through x (lines of the matrix)
            // for (int i = 0; i < N; i++)
            //     thomasAlgorithm(la, lb, lc, c_prime, d_prime, N, RHS[i]);

            // copyMatrices(RHS, V, N);

            // !================================================!
            // ! Calcula V em n + 1/2 -> Resultado vai para RHS !
            // !================================================!
            for (int j = 0; j < N; j++)
            {
                for (int i = 0; i < N; i++)
                {   
                    real x = i * delta_x;
                    real y = j * delta_x;
                    d[i] = phi * V[i][lim(j-1,N)] + (1-2*phi) * V[i][j] + phi * V[i][lim(j+1,N)] + 0.5 * delta_t * forcingTerm(x, y, actualTime);
                }
                
                tridiag(la, lb, lc, c_prime, d_prime, N, d, result);
                for (int i = 0; i < N; i++)
                {
                    RHS[i][j] = result[i];
                }
            }

            // !================================================!
            // ! Calcula V em n + 1 -> Resultado vai para V     !
            // !================================================!
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    real x = i * delta_x;
                    real y = j * delta_x;
                    d[j] = phi * RHS[lim(i-1,N)][j] + (1-2*phi) * RHS[i][j] + phi * RHS[lim(i+1,N)][j] + 0.5 * delta_t * forcingTerm(x, y, actualTime+delta_t);
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

    #ifdef PARALLEL
    // Get (v - solution)²
    errorXerror<<<GRID_SIZE, BLOCK_SIZE>>>(d_V, d_RHS, N, actualTime, delta_x);
    
    // Allocate 2D array for variable
    real *temp = (real *)malloc(N * N * sizeof(real));

    // Copy memory of d_V from device to host of the matrices (2D arrays)
    CUDA_CALL(cudaMemcpy(temp, d_RHS, N * N * sizeof(real), cudaMemcpyDeviceToHost));

    // Calculate the sum of errors²
    real sum = 0.0;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            sum += temp[i*N+j];
    
    // Calculate norm-2 error
    real norm2error = sqrt(delta_x*delta_x*sum);
    #endif // PARALLEL
    #ifdef SERIAL
    real norm2error = calculateNorm2Error(V, exact, N, actualTime, delta_x);
    #endif

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
            #ifdef PARALLEL
            fprintf(fpLast, "%e ", temp[i * N + j]);
            #endif // PARALLEL
            #ifdef SERIAL
            fprintf(fpLast, "%e ", V[i][j]);
            fprintf(fpExact, "%e ", exact[i][j]);
            fprintf(fpErrors, "%e ", abs(V[i][j] - exact[i][j]));
            #endif // SERIAL
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
    #ifdef PARALLEL
    fprintf(fpInfos, "\nFor 1st Part and Transpose -> Grid size %d, Block size %d\n", GRID_SIZE, BLOCK_SIZE);
    fprintf(fpInfos, "Total threads: %d\n", GRID_SIZE*BLOCK_SIZE);
    fprintf(fpInfos, "\nFor 2nd Part -> Grid size: %d, Block size: %d\n", numBlocks, blockSize);
    fprintf(fpInfos, "Total threads: %d\n", numBlocks*blockSize);
    #endif // PARALLEL
    fprintf(fpInfos, "\nNorm-2 Error = %lf\n", norm2error);

    // Close files
    fclose(fpInfos);

    // Free memory
    free(time);

    // Free memory from host
    #ifdef PARALLEL
    free(temp);
    #endif // PARALLEL
    #ifdef SERIAL
    for (int i = 0; i < N; i++)
    {
        free(V[i]);
        free(Rv[i]);
        free(RHS[i]);
        free(exact[i]);
    }
    free(V);
    free(Rv);
    free(RHS);
    free(exact);
    free(c_prime);
    free(d_prime);
    free(d);
    free(result);
    #endif // SERIAL
    free(la);
    free(lb);
    free(lc);
    free(pathToSaveData);
    free(aux);

    #ifdef PARALLEL
    // Free memory from device
    CUDA_CALL(cudaFree(d_V));
    CUDA_CALL(cudaFree(d_RHS));
    CUDA_CALL(cudaFree(d_Rv));
    CUDA_CALL(cudaFree(d_la));
    CUDA_CALL(cudaFree(d_lb));
    CUDA_CALL(cudaFree(d_lc));

    // Reset device
    CUDA_CALL(cudaDeviceReset());
    #endif // PARALLEL

    printf("Simulation finished!\n");

    return;
}


#endif // CONVERGENCE_METHODS_H