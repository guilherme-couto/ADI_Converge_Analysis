#ifndef GPU_METHODS_H
#define GPU_METHODS_H

#include "auxfuncs.h"

void runSimulationGPU(char *method, real delta_t, real delta_x, real theta)
{
    // Number of steps
    int N = round(L / delta_x) + 1;     // Spatial steps (square tissue)
    int M = round(totalTime / delta_t); // Number of time steps

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
    initialize2DVariableWithValue(V, N, V_init);
#endif // CONVERGENCE_ANALYSIS

#ifdef MONODOMAIN
    #ifdef AFHN
    real *W = (real *)malloc(N * N * sizeof(real));
    initialize2DVariableWithValue(W, N, W_init);
    #endif // AFHN
    #ifdef TT2
    real *X_r1 = (real *)malloc(N * N * sizeof(real));
    real *X_r2 = (real *)malloc(N * N * sizeof(real));
    real *X_s = (real *)malloc(N * N * sizeof(real));
    real *m = (real *)malloc(N * N * sizeof(real));
    real *h = (real *)malloc(N * N * sizeof(real));
    real *j = (real *)malloc(N * N * sizeof(real));
    real *d = (real *)malloc(N * N * sizeof(real));
    real *f = (real *)malloc(N * N * sizeof(real));
    real *f2 = (real *)malloc(N * N * sizeof(real));
    real *fCaSS = (real *)malloc(N * N * sizeof(real));
    real *s = (real *)malloc(N * N * sizeof(real));
    real *r = (real *)malloc(N * N * sizeof(real));
    real *Ca_i = (real *)malloc(N * N * sizeof(real));
    real *Ca_SR = (real *)malloc(N * N * sizeof(real));
    real *Ca_SS = (real *)malloc(N * N * sizeof(real));
    real *R_prime = (real *)malloc(N * N * sizeof(real));
    real *Na_i = (real *)malloc(N * N * sizeof(real));
    real *K_i = (real *)malloc(N * N * sizeof(real));
    initialize2DVariableWithValue(X_r1, N, X_r1_init);
    initialize2DVariableWithValue(X_r2, N, X_r2_init);
    initialize2DVariableWithValue(X_s, N, X_s_init);
    initialize2DVariableWithValue(m, N, m_init);
    initialize2DVariableWithValue(h, N, h_init);
    initialize2DVariableWithValue(j, N, j_init);
    initialize2DVariableWithValue(d, N, d_init);
    initialize2DVariableWithValue(f, N, f_init);
    initialize2DVariableWithValue(f2, N, f2_init);
    initialize2DVariableWithValue(fCaSS, N, fCaSS_init);
    initialize2DVariableWithValue(s, N, s_init);
    initialize2DVariableWithValue(r, N, r_init);
    initialize2DVariableWithValue(Ca_i, N, Ca_i_init);
    initialize2DVariableWithValue(Ca_SR, N, Ca_SR_init);
    initialize2DVariableWithValue(Ca_SS, N, Ca_SS_init);
    initialize2DVariableWithValue(R_prime, N, R_prime_init);
    initialize2DVariableWithValue(Na_i, N, Na_i_init);
    initialize2DVariableWithValue(K_i, N, K_i_init);
    #endif // TT2
#endif // MONODOMAIN

#ifdef INIT_WITH_SPIRAL
    char* reference_dt = "0.00010";
    char* reference_dx = "0.00500";
    real real_ref_dx = 0.005f;
    char *pathToSpiralFiles = (char *)malloc(MAX_STRING_SIZE * sizeof(char));
    snprintf(pathToSpiralFiles, MAX_STRING_SIZE * sizeof(char), "./spiral_files/%s/%s/%s/lastV_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize2DVariableFromFile(V, N, pathToSpiralFiles, delta_x, "V", real_ref_dx);
#ifdef AFHN
    snprintf(pathToSpiralFiles, MAX_STRING_SIZE * sizeof(char), "./spiral_files/%s/%s/%s/lastW_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize2DVariableFromFile(W, N, pathToSpiralFiles, delta_x, "W", real_ref_dx);
#endif // AFHN
#ifdef TT2
    snprintf(pathToSpiralFiles, MAX_STRING_SIZE * sizeof(char), "./spiral_files/%s/%s/%s/lastX_r1_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize2DVariableFromFile(X_r1, N, pathToSpiralFiles, delta_x, "X_r1", real_ref_dx);
    snprintf(pathToSpiralFiles, MAX_STRING_SIZE * sizeof(char), "./spiral_files/%s/%s/%s/lastX_r2_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize2DVariableFromFile(X_r2, N, pathToSpiralFiles, delta_x, "X_r2", real_ref_dx);
    snprintf(pathToSpiralFiles, MAX_STRING_SIZE * sizeof(char), "./spiral_files/%s/%s/%s/lastX_s_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize2DVariableFromFile(X_s, N, pathToSpiralFiles, delta_x, "X_s", real_ref_dx);
    snprintf(pathToSpiralFiles, MAX_STRING_SIZE * sizeof(char), "./spiral_files/%s/%s/%s/lastm_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize2DVariableFromFile(m, N, pathToSpiralFiles, delta_x, "m", real_ref_dx);
    snprintf(pathToSpiralFiles, MAX_STRING_SIZE * sizeof(char), "./spiral_files/%s/%s/%s/lasth_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize2DVariableFromFile(h, N, pathToSpiralFiles, delta_x, "h", real_ref_dx);
    snprintf(pathToSpiralFiles, MAX_STRING_SIZE * sizeof(char), "./spiral_files/%s/%s/%s/lastj_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize2DVariableFromFile(j, N, pathToSpiralFiles, delta_x, "j", real_ref_dx);
    snprintf(pathToSpiralFiles, MAX_STRING_SIZE * sizeof(char), "./spiral_files/%s/%s/%s/lastd_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize2DVariableFromFile(d, N, pathToSpiralFiles, delta_x, "d", real_ref_dx);
    snprintf(pathToSpiralFiles, MAX_STRING_SIZE * sizeof(char), "./spiral_files/%s/%s/%s/lastf_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize2DVariableFromFile(f, N, pathToSpiralFiles, delta_x, "f", real_ref_dx);
    snprintf(pathToSpiralFiles, MAX_STRING_SIZE * sizeof(char), "./spiral_files/%s/%s/%s/lastf2_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize2DVariableFromFile(f2, N, pathToSpiralFiles, delta_x, "f2", real_ref_dx);
    snprintf(pathToSpiralFiles, MAX_STRING_SIZE * sizeof(char), "./spiral_files/%s/%s/%s/lastfCaSS_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize2DVariableFromFile(fCaSS, N, pathToSpiralFiles, delta_x, "fCaSS", real_ref_dx);
    snprintf(pathToSpiralFiles, MAX_STRING_SIZE * sizeof(char), "./spiral_files/%s/%s/%s/lasts_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize2DVariableFromFile(s, N, pathToSpiralFiles, delta_x, "s", real_ref_dx);
    snprintf(pathToSpiralFiles, MAX_STRING_SIZE * sizeof(char), "./spiral_files/%s/%s/%s/lastr_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize2DVariableFromFile(r, N, pathToSpiralFiles, delta_x, "r", real_ref_dx);
    snprintf(pathToSpiralFiles, MAX_STRING_SIZE * sizeof(char), "./spiral_files/%s/%s/%s/lastCa_i_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize2DVariableFromFile(Ca_i, N, pathToSpiralFiles, delta_x, "Ca_i", real_ref_dx);
    snprintf(pathToSpiralFiles, MAX_STRING_SIZE * sizeof(char), "./spiral_files/%s/%s/%s/lastCa_SR_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize2DVariableFromFile(Ca_SR, N, pathToSpiralFiles, delta_x, "Ca_SR", real_ref_dx);
    snprintf(pathToSpiralFiles, MAX_STRING_SIZE * sizeof(char), "./spiral_files/%s/%s/%s/lastCa_SS_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize2DVariableFromFile(Ca_SS, N, pathToSpiralFiles, delta_x, "Ca_SS", real_ref_dx);
    snprintf(pathToSpiralFiles, MAX_STRING_SIZE * sizeof(char), "./spiral_files/%s/%s/%s/lastR_prime_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize2DVariableFromFile(R_prime, N, pathToSpiralFiles, delta_x, "R_prime", real_ref_dx);
    snprintf(pathToSpiralFiles, MAX_STRING_SIZE * sizeof(char), "./spiral_files/%s/%s/%s/lastNa_i_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize2DVariableFromFile(Na_i, N, pathToSpiralFiles, delta_x, "Na_i", real_ref_dx);
    snprintf(pathToSpiralFiles, MAX_STRING_SIZE * sizeof(char), "./spiral_files/%s/%s/%s/lastK_i_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize2DVariableFromFile(K_i, N, pathToSpiralFiles, delta_x, "K_i", real_ref_dx);
#endif // TT2
    free(pathToSpiralFiles);
#endif // INIT_WITH_SPIRAL

    // Auxiliary arrays for Thomas algorithm
    real *la = (real *)malloc(N * sizeof(real)); // subdiagonal
    real *lb = (real *)malloc(N * sizeof(real)); // diagonal
    real *lc = (real *)malloc(N * sizeof(real)); // superdiagonal

    // Populate auxiliary arrays for Thomas algorithm
    real phi = delta_t / (delta_x * delta_x);
    if (strcmp(method, "ADI") == 0 || strcmp(method, "SSI-ADI") == 0)
    {
#ifdef AFHN
        populateDiagonalThomasAlgorithm(la, lb, lc, N, 0.5f * phi * (sigma / (Cm * chi)));
#endif // AFHN
#ifdef TT2
        populateDiagonalThomasAlgorithm(la, lb, lc, N, 0.5f * phi * (sigma / (chi)));
#endif // TT2
    }
    else if (strstr(method, "theta") != NULL)
    {
#ifdef AFHN
        populateDiagonalThomasAlgorithm(la, lb, lc, N, theta * phi * (sigma / (Cm * chi)));
#endif // AFHN
#ifdef TT2
        populateDiagonalThomasAlgorithm(la, lb, lc, N, theta * phi * (sigma / chi));
#endif // TT2
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

#ifdef MONODOMAIN
#ifdef AFHN
    real *d_W;
    CUDA_CALL(cudaMalloc(&d_W, N * N * sizeof(real)));
    CUDA_CALL(cudaMemcpy(d_W, W, N * N * sizeof(real), cudaMemcpyHostToDevice));
#endif // AFHN
#ifdef TT2
    real *d_X_r1, *d_X_r2, *d_X_s, *d_m, *d_h, *d_j, *d_d, *d_f, *d_f2, *d_fCaSS, *d_s, *d_r, *d_Ca_i, *d_Ca_SR, *d_Ca_SS, *d_R_prime, *d_Na_i, *d_K_i;
    CUDA_CALL(cudaMalloc(&d_X_r1, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_X_r2, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_X_s, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_m, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_h, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_j, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_d, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_f, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_f2, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_fCaSS, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_s, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_r, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_Ca_i, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_Ca_SR, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_Ca_SS, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_R_prime, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_Na_i, N * N * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_K_i, N * N * sizeof(real)));
    CUDA_CALL(cudaMemcpy(d_X_r1, X_r1, N * N * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_X_r2, X_r2, N * N * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_X_s, X_s, N * N * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_m, m, N * N * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_h, h, N * N * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_j, j, N * N * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_d, d, N * N * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_f, f, N * N * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_f2, f2, N * N * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_fCaSS, fCaSS, N * N * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_s, s, N * N * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_r, r, N * N * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_Ca_i, Ca_i, N * N * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_Ca_SR, Ca_SR, N * N * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_Ca_SS, Ca_SS, N * N * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_R_prime, R_prime, N * N * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_Na_i, Na_i, N * N * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_K_i, K_i, N * N * sizeof(real), cudaMemcpyHostToDevice));
#endif // TT2
#endif // MONODOMAIN

#ifndef CONVERGENCE_ANALYSIS
#ifdef MONODOMAIN
    // Allocate array for the stimuli
    Stimulus *stimuli = (Stimulus *)malloc(numberOfStimuli * sizeof(Stimulus));
    populateStimuli(stimuli, delta_x);

    // Allocate and copy array for the stimuli on device
    Stimulus *d_stimuli;
    CUDA_CALL(cudaMalloc(&d_stimuli, numberOfStimuli * sizeof(Stimulus)));
    CUDA_CALL(cudaMemcpy(d_stimuli, stimuli, numberOfStimuli * sizeof(Stimulus), cudaMemcpyHostToDevice));
#endif // MONODOMAIN
#endif // not CONVERGENCE_ANALYSIS

    // CUDA grid and block allocation
    // For Thomas algorithm kernel
    int numBlocks = N / 100;
    if (numBlocks == 0)
        numBlocks = 1;

    int blockSize = round(N / numBlocks) + 1;
    if (blockSize % 32 != 0)
        blockSize = 32 * ((blockSize / 32) + 1);

    // For ODEs and Transpose kernels
    int GRID_SIZE = ceil((N * N * 1.0f) / (BLOCK_SIZE * 1.0f));

    // Create directories
    char *pathToSaveData = (char *)malloc(MAX_STRING_SIZE * sizeof(char));
    createDirectories(method, theta, pathToSaveData);

#ifdef SAVE_FRAMES
    // Save frames
    char framesPath[MAX_STRING_SIZE];
    FILE *fpFrames;
    snprintf(framesPath, MAX_STRING_SIZE * sizeof(char), "%s/frames/frames_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);
    fpFrames = fopen(framesPath, "w");
#endif // SAVE_FRAMES

    // Measure velocity
    real stim_velocity = 0.0;
    real first_point_time = 0.0;
    real last_point_time = 0.0;
    bool aux_stim_velocity_flag = false;
    bool stim_velocity_measured = false;

    int timeStepCounter = 0;
    real actualTime = 0.0f;

    // Measure total execution time
    real startTime = omp_get_wtime();

    if (strcmp(method, "ADI") == 0)
    {
        while (timeStepCounter < M)
        {
            // Get time step
            actualTime = time[timeStepCounter];

            printf("ADI for GPU not implemented yet\n");
            exit(EXIT_FAILURE);

            // Update time step counter
            timeStepCounter++;
        }
    }

#if defined(LINMONO) || defined(MONODOMAIN)
    else if (strcmp(method, "SSI-ADI") == 0 || strcmp(method, "theta-ADI") == 0)
    {
        while (timeStepCounter < M)
        {
            // Get time step
            actualTime = time[timeStepCounter];

            // ================================================!
            //  Calculate Approx. and ODEs                     !
            // ================================================!

#if defined(CONVERGENCE_ANALYSIS) && defined(AFHN)
            computeApproxSSI<<<GRID_SIZE, BLOCK_SIZE>>>(N, delta_t, phi, delta_x, actualTime, d_V, d_Vtilde, d_partRHS, d_W);
#elif defined(AFHN)
            computeApprox<<<GRID_SIZE, BLOCK_SIZE>>>(N, delta_t, phi, delta_x, actualTime, d_V, d_partRHS, d_W, d_stimuli);
#endif // CONVERGENCE_ANALYSIS
#ifdef TT2
            computeApprox<<<GRID_SIZE, BLOCK_SIZE>>>(N, delta_t, phi, delta_x, actualTime, d_V, d_partRHS, d_X_r1, d_X_r2, d_X_s, d_m, d_h, d_j, d_d, d_f, d_f2, d_fCaSS, d_s, d_r, d_Ca_i, d_Ca_SR, d_Ca_SS, d_R_prime, d_Na_i, d_K_i, d_stimuli);
#endif // TT2
            cudaDeviceSynchronize();

            // ================================================!
            //  Calculate V at n+1/2 -> Result goes to RHS     !
            //  diffusion implicit in y and explicit in x      !
            // ================================================!
            (strcmp(method, "SSI-ADI") == 0)
                ? prepareRHSwithjDiff<<<GRID_SIZE, BLOCK_SIZE>>>(N, phi, 0.5f, d_V, d_RHS, d_partRHS)
                : prepareRHSwithjDiff<<<GRID_SIZE, BLOCK_SIZE>>>(N, phi, (1.0f - theta), d_V, d_RHS, d_partRHS);
            cudaDeviceSynchronize();

            parallelThomas<<<numBlocks, blockSize>>>(N, d_RHS, d_la, d_lb, d_lc);
            cudaDeviceSynchronize();

            // =========================================================!
            //  Transpose d_RHS to d_Vtilde (d_Vtilde as an aux var)    !
            // =========================================================!
            parallelTranspose<<<GRID_SIZE, BLOCK_SIZE>>>(N, d_RHS, d_Vtilde);
            cudaDeviceSynchronize();

            // ================================================!
            //  Calculate V at n+1 -> Result goes to V         !
            //  diffusion implicit in x and explicit in y      !
            // ================================================!
            (strcmp(method, "SSI-ADI") == 0)
                ? prepareRHSwithiDiff<<<GRID_SIZE, BLOCK_SIZE>>>(N, phi, 0.5f, d_Vtilde, d_V, d_partRHS)
                : prepareRHSwithiDiff<<<GRID_SIZE, BLOCK_SIZE>>>(N, phi, (1.0f - theta), d_Vtilde, d_V, d_partRHS);
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

            // Calculate stim velocity
            if (!stim_velocity_measured)
            {
                // Copy memory of d_V from device to host V
                CUDA_CALL(cudaMemcpy(V, d_V, N * N * sizeof(real), cudaMemcpyDeviceToHost));
                cudaDeviceSynchronize();

                real begin = L / 3.0f;
                real end = 2.0f * begin;

                if (!aux_stim_velocity_flag)
                {
                    int first_point_index = round(begin / delta_x) + 1;
                    if (V[first_point_index] > 10.0)
                    {
                        first_point_time = actualTime;
                        aux_stim_velocity_flag = true;
                    }
                }
                else
                {
                    int last_point_index = round(end / delta_x) + 1;
                    if (V[last_point_index] > 10.0)
                    {
                        last_point_time = actualTime;
                        stim_velocity = (end - begin) / (last_point_time - first_point_time); // cm/ms
                        stim_velocity = stim_velocity * 10.0; // m/s
                        stim_velocity_measured = true;
                        printf("Stim velocity (measured from %f to %f cm) is %lf m/s\n", begin, end, stim_velocity);
                    }
                }
            }

            // Update time step counter
            timeStepCounter++;
        }
    }
#endif // LINMONO || MONODOMAIN

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
    real norm2error = calculateNorm2Error(V, exact, N, totalTime, delta_x);
#endif // CONVERGENCE_ANALYSIS

    // Write infos to file
    char infosFilePath[MAX_STRING_SIZE];
    snprintf(infosFilePath, MAX_STRING_SIZE * sizeof(char), "%s/infos/infos_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);
    FILE *fpInfos = fopen(infosFilePath, "w");
    printf("Infos saved to %s\n", infosFilePath);
    fprintf(fpInfos, "Domain Length = %d, Time = %f\n", L, totalTime);
    fprintf(fpInfos, "delta_x = %lf, Space steps N = %d, N*N = %d\n", delta_x, N, N * N);
    fprintf(fpInfos, "delta_t = %lf, Time steps = %d\n", delta_t, M);
    fprintf(fpInfos, "Stimulus velocity = %lf m/s\n", stim_velocity);
    fprintf(fpInfos, "Method %s\n", method);
    fprintf(fpInfos, "\nFor Approx, ODEs and Transpose -> Grid size %d, Block size %d\n", GRID_SIZE, BLOCK_SIZE);
    fprintf(fpInfos, "Total threads: %d\n", GRID_SIZE * BLOCK_SIZE);
    fprintf(fpInfos, "\nFor Thomas -> Grid size: %d, Block size: %d\n", numBlocks, blockSize);
    fprintf(fpInfos, "Total threads: %d\n", numBlocks * blockSize);
    fprintf(fpInfos, "\nTotal execution time = %lf\n", elapsedTime);
#ifdef CONVERGENCE_ANALYSIS
    fprintf(fpInfos, "\nNorm-2 Error = %lf\n", norm2error);
#endif // CONVERGENCE_ANALYSIS
    fclose(fpInfos);

    // Save last frame
    char lastFrameFilePath[MAX_STRING_SIZE];
    snprintf(lastFrameFilePath, MAX_STRING_SIZE * sizeof(char), "%s/lastframe/last_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);
    FILE *fpLast = fopen(lastFrameFilePath, "w");
    printf("Last frame saved to %s\n", lastFrameFilePath);

#ifdef SAVE_LAST_STATE
#ifdef AFHN
    char lastFrameFilePathV[MAX_STRING_SIZE], lastFrameFilePathW[MAX_STRING_SIZE];

    snprintf(lastFrameFilePathV, MAX_STRING_SIZE * sizeof(char), "%s/lastframe/lastV_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);
    snprintf(lastFrameFilePathW, MAX_STRING_SIZE * sizeof(char), "%s/lastframe/lastW_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);

    FILE *fpLastV = fopen(lastFrameFilePathV, "w");
    printf("Last V frame saved to %s\n", lastFrameFilePathV);
    FILE *fpLastW = fopen(lastFrameFilePathW, "w");
    printf("Last W frame saved to %s\n", lastFrameFilePathW);
    CUDA_CALL(cudaMemcpy(W, d_W, N * N * sizeof(real), cudaMemcpyDeviceToHost));
#endif // AFHN
#ifdef TT2
    char lastFrameFilePathV[MAX_STRING_SIZE], lastFrameFilePathX_r1[MAX_STRING_SIZE], lastFrameFilePathX_r2[MAX_STRING_SIZE], lastFrameFilePathX_s[MAX_STRING_SIZE], lastFrameFilePathm[MAX_STRING_SIZE], lastFrameFilePathh[MAX_STRING_SIZE], lastFrameFilePathj[MAX_STRING_SIZE], lastFrameFilePathd[MAX_STRING_SIZE], lastFrameFilePathf[MAX_STRING_SIZE], lastFrameFilePathf2[MAX_STRING_SIZE], lastFrameFilePathfCaSS[MAX_STRING_SIZE], lastFrameFilePaths[MAX_STRING_SIZE], lastFrameFilePathr[MAX_STRING_SIZE], lastFrameFilePathR_prime[MAX_STRING_SIZE], lastFrameFilePathCa_i[MAX_STRING_SIZE], lastFrameFilePathCa_SR[MAX_STRING_SIZE], lastFrameFilePathCa_SS[MAX_STRING_SIZE], lastFrameFilePathNa_i[MAX_STRING_SIZE], lastFrameFilePathK_i[MAX_STRING_SIZE];
    
    snprintf(lastFrameFilePathV, MAX_STRING_SIZE * sizeof(char), "%s/lastframe/lastV_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);
    snprintf(lastFrameFilePathX_r1, MAX_STRING_SIZE * sizeof(char), "%s/lastframe/lastX_r1_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);
    snprintf(lastFrameFilePathX_r2, MAX_STRING_SIZE * sizeof(char), "%s/lastframe/lastX_r2_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);
    snprintf(lastFrameFilePathX_s, MAX_STRING_SIZE * sizeof(char), "%s/lastframe/lastX_s_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);
    snprintf(lastFrameFilePathm, MAX_STRING_SIZE * sizeof(char), "%s/lastframe/lastm_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);
    snprintf(lastFrameFilePathh, MAX_STRING_SIZE * sizeof(char), "%s/lastframe/lasth_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);
    snprintf(lastFrameFilePathj, MAX_STRING_SIZE * sizeof(char), "%s/lastframe/lastj_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);
    snprintf(lastFrameFilePathd, MAX_STRING_SIZE * sizeof(char), "%s/lastframe/lastd_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);
    snprintf(lastFrameFilePathf, MAX_STRING_SIZE * sizeof(char), "%s/lastframe/lastf_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);
    snprintf(lastFrameFilePathf2, MAX_STRING_SIZE * sizeof(char), "%s/lastframe/lastf2_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);
    snprintf(lastFrameFilePathfCaSS, MAX_STRING_SIZE * sizeof(char), "%s/lastframe/lastfCaSS_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);
    snprintf(lastFrameFilePaths, MAX_STRING_SIZE * sizeof(char), "%s/lastframe/lasts_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);
    snprintf(lastFrameFilePathr, MAX_STRING_SIZE * sizeof(char), "%s/lastframe/lastr_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);
    snprintf(lastFrameFilePathR_prime, MAX_STRING_SIZE * sizeof(char), "%s/lastframe/lastR_prime_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);
    snprintf(lastFrameFilePathCa_i, MAX_STRING_SIZE * sizeof(char), "%s/lastframe/lastCa_i_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);
    snprintf(lastFrameFilePathCa_SR, MAX_STRING_SIZE * sizeof(char), "%s/lastframe/lastCa_SR_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);
    snprintf(lastFrameFilePathCa_SS, MAX_STRING_SIZE * sizeof(char), "%s/lastframe/lastCa_SS_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);
    snprintf(lastFrameFilePathNa_i, MAX_STRING_SIZE * sizeof(char), "%s/lastframe/lastNa_i_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);
    snprintf(lastFrameFilePathK_i, MAX_STRING_SIZE * sizeof(char), "%s/lastframe/lastK_i_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);

    FILE *fpLastV = fopen(lastFrameFilePathV, "w");
    printf("Last V frame saved to %s\n", lastFrameFilePathV);
    FILE *fpLastX_r1 = fopen(lastFrameFilePathX_r1, "w");
    printf("Last X_r1 frame saved to %s\n", lastFrameFilePathX_r1);
    CUDA_CALL(cudaMemcpy(X_r1, d_X_r1, N * N * sizeof(real), cudaMemcpyDeviceToHost));
    FILE *fpLastX_r2 = fopen(lastFrameFilePathX_r2, "w");
    printf("Last X_r2 frame saved to %s\n", lastFrameFilePathX_r2);
    CUDA_CALL(cudaMemcpy(X_r2, d_X_r2, N * N * sizeof(real), cudaMemcpyDeviceToHost));
    FILE *fpLastX_s = fopen(lastFrameFilePathX_s, "w");
    printf("Last X_s frame saved to %s\n", lastFrameFilePathX_s);
    CUDA_CALL(cudaMemcpy(X_s, d_X_s, N * N * sizeof(real), cudaMemcpyDeviceToHost));
    FILE *fpLastm = fopen(lastFrameFilePathm, "w");
    printf("Last m frame saved to %s\n", lastFrameFilePathm);
    CUDA_CALL(cudaMemcpy(m, d_m, N * N * sizeof(real), cudaMemcpyDeviceToHost));
    FILE *fpLasth = fopen(lastFrameFilePathh, "w");
    printf("Last h frame saved to %s\n", lastFrameFilePathh);
    CUDA_CALL(cudaMemcpy(h, d_h, N * N * sizeof(real), cudaMemcpyDeviceToHost));
    FILE *fpLastj = fopen(lastFrameFilePathj, "w");
    printf("Last j frame saved to %s\n", lastFrameFilePathj);
    CUDA_CALL(cudaMemcpy(j, d_j, N * N * sizeof(real), cudaMemcpyDeviceToHost));
    FILE *fpLastd = fopen(lastFrameFilePathd, "w");
    printf("Last d frame saved to %s\n", lastFrameFilePathd);
    CUDA_CALL(cudaMemcpy(d, d_d, N * N * sizeof(real), cudaMemcpyDeviceToHost));
    FILE *fpLastf = fopen(lastFrameFilePathf, "w");
    printf("Last f frame saved to %s\n", lastFrameFilePathf);
    CUDA_CALL(cudaMemcpy(f, d_f, N * N * sizeof(real), cudaMemcpyDeviceToHost));
    FILE *fpLastf2 = fopen(lastFrameFilePathf2, "w");
    printf("Last f2 frame saved to %s\n", lastFrameFilePathf2);
    CUDA_CALL(cudaMemcpy(f2, d_f2, N * N * sizeof(real), cudaMemcpyDeviceToHost));
    FILE *fpLastfCaSS = fopen(lastFrameFilePathfCaSS, "w");
    printf("Last fCaSS frame saved to %s\n", lastFrameFilePathfCaSS);
    CUDA_CALL(cudaMemcpy(fCaSS, d_fCaSS, N * N * sizeof(real), cudaMemcpyDeviceToHost));
    FILE *fpLasts = fopen(lastFrameFilePaths, "w");
    printf("Last s frame saved to %s\n", lastFrameFilePaths);
    CUDA_CALL(cudaMemcpy(s, d_s, N * N * sizeof(real), cudaMemcpyDeviceToHost));
    FILE *fpLastr = fopen(lastFrameFilePathr, "w");
    printf("Last r frame saved to %s\n", lastFrameFilePathr);
    CUDA_CALL(cudaMemcpy(r, d_r, N * N * sizeof(real), cudaMemcpyDeviceToHost));
    FILE *fpLastR_prime = fopen(lastFrameFilePathR_prime, "w");
    printf("Last R_prime frame saved to %s\n", lastFrameFilePathR_prime);
    CUDA_CALL(cudaMemcpy(R_prime, d_R_prime, N * N * sizeof(real), cudaMemcpyDeviceToHost));
    FILE *fpLastCa_i = fopen(lastFrameFilePathCa_i, "w");
    printf("Last Ca_i frame saved to %s\n", lastFrameFilePathCa_i);
    CUDA_CALL(cudaMemcpy(Ca_i, d_Ca_i, N * N * sizeof(real), cudaMemcpyDeviceToHost));
    FILE *fpLastCa_SR = fopen(lastFrameFilePathCa_SR, "w");
    printf("Last Ca_SR frame saved to %s\n", lastFrameFilePathCa_SR);
    CUDA_CALL(cudaMemcpy(Ca_SR, d_Ca_SR, N * N * sizeof(real), cudaMemcpyDeviceToHost));
    FILE *fpLastCa_SS = fopen(lastFrameFilePathCa_SS, "w");
    printf("Last Ca_SS frame saved to %s\n", lastFrameFilePathCa_SS);
    CUDA_CALL(cudaMemcpy(Ca_SS, d_Ca_SS, N * N * sizeof(real), cudaMemcpyDeviceToHost));
    FILE *fpLastNa_i = fopen(lastFrameFilePathNa_i, "w");
    printf("Last Na_i frame saved to %s\n", lastFrameFilePathNa_i);
    CUDA_CALL(cudaMemcpy(Na_i, d_Na_i, N * N * sizeof(real), cudaMemcpyDeviceToHost));
    FILE *fpLastK_i = fopen(lastFrameFilePathK_i, "w");
    printf("Last K_i frame saved to %s\n", lastFrameFilePathK_i);
    CUDA_CALL(cudaMemcpy(K_i, d_K_i, N * N * sizeof(real), cudaMemcpyDeviceToHost));
#endif // TT2
#endif // SAVE_LAST_STATE
             
#ifdef CONVERGENCE_ANALYSIS
    char exactFilePath[MAX_STRING_SIZE];
    snprintf(exactFilePath, MAX_STRING_SIZE * sizeof(char), "%s/exact/exact_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);
    FILE *fpExact = fopen(exactFilePath, "w");
    printf("Exact solution saved to %s\n", exactFilePath);
    char errorsFilePath[MAX_STRING_SIZE];
    snprintf(errorsFilePath, MAX_STRING_SIZE * sizeof(char), "%s/errors/errors_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);
    FILE *fpErrors = fopen(errorsFilePath, "w");
    printf("Errors saved to %s\n", errorsFilePath);
#endif // CONVERGENCE_ANALYSIS
    for (int i = 0; i < N; i++)
    {
        for (int _j = 0; _j < N; _j++)
        {
            int index = i * N + _j;
            fprintf(fpLast, "%e ", V[index]);
#ifdef CONVERGENCE_ANALYSIS
            fprintf(fpExact, "%e ", exact[i][_j]);
            fprintf(fpErrors, "%e ", abs(V[index] - exact[i][_j]));
#endif // CONVERGENCE_ANALYSIS
#ifdef SAVE_LAST_STATE
#ifdef AFHN
            fprintf(fpLastV, "%e ", V[index]);
            fprintf(fpLastW, "%e ", W[index]);
#endif // AFHN
#ifdef TT2
            fprintf(fpLastV, "%e ", V[index]);
            fprintf(fpLastX_r1, "%e ", X_r1[index]);
            fprintf(fpLastX_r2, "%e ", X_r2[index]);
            fprintf(fpLastX_s, "%e ", X_s[index]);
            fprintf(fpLastm, "%e ", m[index]);
            fprintf(fpLasth, "%e ", h[index]);
            fprintf(fpLastj, "%e ", j[index]);
            fprintf(fpLastd, "%e ", d[index]);
            fprintf(fpLastf, "%e ", f[index]);
            fprintf(fpLastf2, "%e ", f2[index]);
            fprintf(fpLastfCaSS, "%e ", fCaSS[index]);
            fprintf(fpLasts, "%e ", s[index]);
            fprintf(fpLastr, "%e ", r[index]);
            fprintf(fpLastR_prime, "%e ", R_prime[index]);
            fprintf(fpLastCa_i, "%e ", Ca_i[index]);
            fprintf(fpLastCa_SR, "%e ", Ca_SR[index]);
            fprintf(fpLastCa_SS, "%e ", Ca_SS[index]);
            fprintf(fpLastNa_i, "%e ", Na_i[index]);
            fprintf(fpLastK_i, "%e ", K_i[index]);
#endif // TT2
#endif // SAVE_LAST_STATE
        }
        fprintf(fpLast, "\n");
#ifdef CONVERGENCE_ANALYSIS
        fprintf(fpExact, "\n");
        fprintf(fpErrors, "\n");
#endif // CONVERGENCE_ANALYSIS
#ifdef SAVE_LAST_STATE
        #ifdef AFHN
        fprintf(fpLastV, "\n");
        fprintf(fpLastW, "\n");
        #endif // AFHN
        #ifdef TT2
        fprintf(fpLastV, "\n");
        fprintf(fpLastX_r1, "\n");
        fprintf(fpLastX_r2, "\n");
        fprintf(fpLastX_s, "\n");
        fprintf(fpLastm, "\n");
        fprintf(fpLasth, "\n");
        fprintf(fpLastj, "\n");
        fprintf(fpLastd, "\n");
        fprintf(fpLastf, "\n");
        fprintf(fpLastf2, "\n");
        fprintf(fpLastfCaSS, "\n");
        fprintf(fpLasts, "\n");
        fprintf(fpLastr, "\n");
        fprintf(fpLastR_prime, "\n");
        fprintf(fpLastCa_i, "\n");
        fprintf(fpLastCa_SR, "\n");
        fprintf(fpLastCa_SS, "\n");
        fprintf(fpLastNa_i, "\n");
        fprintf(fpLastK_i, "\n");
        #endif // TT2
#endif // SAVE_LAST_STATE
    }
    fclose(fpLast);
#ifdef CONVERGENCE_ANALYSIS
    fclose(fpExact);
    fclose(fpErrors);
#endif // CONVERGENCE_ANALYSIS
#ifdef SAVE_LAST_STATE
    #ifdef AFHN
    fclose(fpLastV);
    fclose(fpLastW);
    #endif // AFHN
    #ifdef TT2
    fclose(fpLastV);
    fclose(fpLastX_r1);
    fclose(fpLastX_r2);
    fclose(fpLastX_s);
    fclose(fpLastm);
    fclose(fpLasth);
    fclose(fpLastj);
    fclose(fpLastd);
    fclose(fpLastf);
    fclose(fpLastf2);
    fclose(fpLastfCaSS);
    fclose(fpLasts);
    fclose(fpLastr);
    fclose(fpLastR_prime);
    fclose(fpLastCa_i);
    fclose(fpLastCa_SR);
    fclose(fpLastCa_SS);
    fclose(fpLastNa_i);
    fclose(fpLastK_i);
    #endif // TT2
#endif // SAVE_LAST_STATE

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
#ifdef MONODOMAIN
    #ifdef AFHN
    free(W);
    #endif // AFHN
    #ifdef TT2
    free(X_r1);
    free(X_r2);
    free(X_s);
    free(m);
    free(h);
    free(j);
    free(d);
    free(f);
    free(f2);
    free(fCaSS);
    free(s);
    free(r);
    free(Ca_i);
    free(Ca_SR);
    free(Ca_SS);
    free(R_prime);
    free(Na_i);
    free(K_i);
    #endif // TT2
    #ifndef CONVERGENCE_ANALYSIS
    free(stimuli);
    #endif // not CONVERGENCE_ANALYSIS
#endif // MONODOMAIN

    // Free memory from device
    CUDA_CALL(cudaFree(d_V));
    CUDA_CALL(cudaFree(d_Vtilde));
    CUDA_CALL(cudaFree(d_RHS));
    CUDA_CALL(cudaFree(d_partRHS));
    CUDA_CALL(cudaFree(d_la));
    CUDA_CALL(cudaFree(d_lb));
    CUDA_CALL(cudaFree(d_lc));
#ifdef MONODOMAIN
#ifdef AFHN
    CUDA_CALL(cudaFree(d_W));
#endif // AFHN
#ifdef TT2
    CUDA_CALL(cudaFree(d_X_r1));
    CUDA_CALL(cudaFree(d_X_r2));
    CUDA_CALL(cudaFree(d_X_s));
    CUDA_CALL(cudaFree(d_m));
    CUDA_CALL(cudaFree(d_h));
    CUDA_CALL(cudaFree(d_j));
    CUDA_CALL(cudaFree(d_d));
    CUDA_CALL(cudaFree(d_f));
    CUDA_CALL(cudaFree(d_f2));
    CUDA_CALL(cudaFree(d_fCaSS));
    CUDA_CALL(cudaFree(d_s));
    CUDA_CALL(cudaFree(d_r));
    CUDA_CALL(cudaFree(d_Ca_i));
    CUDA_CALL(cudaFree(d_Ca_SR));
    CUDA_CALL(cudaFree(d_Ca_SS));
    CUDA_CALL(cudaFree(d_R_prime));
    CUDA_CALL(cudaFree(d_Na_i));
    CUDA_CALL(cudaFree(d_K_i));
#endif // TT2
#ifndef CONVERGENCE_ANALYSIS
    CUDA_CALL(cudaFree(d_stimuli));
#endif // not CONVERGENCE_ANALYSIS
#endif // MONODOMAIN

    // Reset device
    CUDA_CALL(cudaDeviceReset());

    return;
}

#endif // GPU_METHODS_H