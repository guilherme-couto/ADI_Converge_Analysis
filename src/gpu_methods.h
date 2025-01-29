#ifndef GPU_METHODS_H
#define GPU_METHODS_H

#include "auxfuncs.h"

void runSimulationGPU(real delta_t, real delta_x, real delta_y)
{
    // Number of steps
    int M = round(totalTime / delta_t); // Number of time steps
    int Nx = round(Lx / delta_x) + 1;   // Spatial steps in x
    int Ny = round(Ly / delta_y) + 1; // Spatial steps in y

    printf("Nx = %d\n", Nx);
    printf("Ny = %d\n", Ny);
    printf("Points in the domain = %d\n", Nx * Ny);

    // Allocate and populate time array
    real *time = (real *)malloc(M * sizeof(real));
    initializeTimeArray(time, M, delta_t);

    // Allocate arrays for variables (2D matrices will be flattened)
    real *V, *Vtilde, *RHS, *partRHS;
    V = (real *)malloc(Nx * Ny * sizeof(real));
    Vtilde = (real *)malloc(Nx * Ny * sizeof(real));
    RHS = (real *)malloc(Nx * Ny * sizeof(real));
    partRHS = (real *)malloc(Nx * Ny * sizeof(real));

#ifdef AFHN

    real *W = (real *)malloc(Nx * Ny * sizeof(real));

#endif // AFHN
#ifdef TT2

    real *X_r1 = (real *)malloc(Nx * Ny * sizeof(real));
    real *X_r2 = (real *)malloc(Nx * Ny * sizeof(real));
    real *X_s = (real *)malloc(Nx * Ny * sizeof(real));
    real *m = (real *)malloc(Nx * Ny * sizeof(real));
    real *h = (real *)malloc(Nx * Ny * sizeof(real));
    real *_j = (real *)malloc(Nx * Ny * sizeof(real));
    real *d = (real *)malloc(Nx * Ny * sizeof(real));
    real *f = (real *)malloc(Nx * Ny * sizeof(real));
    real *f2 = (real *)malloc(Nx * Ny * sizeof(real));
    real *fCaSS = (real *)malloc(Nx * Ny * sizeof(real));
    real *s = (real *)malloc(Nx * Ny * sizeof(real));
    real *r = (real *)malloc(Nx * Ny * sizeof(real));
    real *Ca_i = (real *)malloc(Nx * Ny * sizeof(real));
    real *Ca_SR = (real *)malloc(Nx * Ny * sizeof(real));
    real *Ca_SS = (real *)malloc(Nx * Ny * sizeof(real));
    real *R_prime = (real *)malloc(Nx * Ny * sizeof(real));
    real *Na_i = (real *)malloc(Nx * Ny * sizeof(real));
    real *K_i = (real *)malloc(Nx * Ny * sizeof(real));

#endif // TT2

#ifdef CONVERGENCE_ANALYSIS_FORCING_TERM

    initialize2DVariableWithExactSolution(V, Nx, Ny, delta_x, delta_y);

#else // not CONVERGENCE_ANALYSIS_FORCING_TERM

    initialize2DVariableWithValue(V, Nx, Ny, V_init);

#endif // CONVERGENCE_ANALYSIS_FORCING_TERM

#ifdef AFHN

    initialize2DVariableWithValue(W, Nx, Ny, W_init);

#endif // AFHN

#ifdef TT2

    initialize2DVariableWithValue(X_r1, Nx, Ny, X_r1_init);
    initialize2DVariableWithValue(X_r2, Nx, Ny, X_r2_init);
    initialize2DVariableWithValue(X_s, Nx, Ny, X_s_init);
    initialize2DVariableWithValue(m, Nx, Ny, m_init);
    initialize2DVariableWithValue(h, Nx, Ny, h_init);
    initialize2DVariableWithValue(_j, Nx, Ny, j_init);
    initialize2DVariableWithValue(d, Nx, Ny, d_init);
    initialize2DVariableWithValue(f, Nx, Ny, f_init);
    initialize2DVariableWithValue(f2, Nx, Ny, f2_init);
    initialize2DVariableWithValue(fCaSS, Nx, Ny, fCaSS_init);
    initialize2DVariableWithValue(s, Nx, Ny, s_init);
    initialize2DVariableWithValue(r, Nx, Ny, r_init);
    initialize2DVariableWithValue(Ca_i, Nx, Ny, Ca_i_init);
    initialize2DVariableWithValue(Ca_SR, Nx, Ny, Ca_SR_init);
    initialize2DVariableWithValue(Ca_SS, Nx, Ny, Ca_SS_init);
    initialize2DVariableWithValue(R_prime, Nx, Ny, R_prime_init);
    initialize2DVariableWithValue(Na_i, Nx, Ny, Na_i_init);
    initialize2DVariableWithValue(K_i, Nx, Ny, K_i_init);

#endif // TT2

#ifdef RESTORE_STATE

    printf("\n");
    printf("Restoring state variables...\n");

    // Initialize variables with a solution
    real real_ref_dx = 0.0005f;
    real real_def_dy = 0.0005f;

    char *pathToRestoreStateFiles = (char *)malloc(MAX_STRING_SIZE * sizeof(char));

    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastframeV.txt", REAL_TYPE, PROBLEM, CELL_MODEL);
    initialize2DVariableFromFile(V, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "V", real_ref_dx, real_def_dy);

#ifdef AFHN

    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastframeW.txt", REAL_TYPE, PROBLEM, CELL_MODEL);
    initialize2DVariableFromFile(W, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "W", real_ref_dx, real_def_dy);

#endif // AFHN

#ifdef TT2
    
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastframeX_r1.txt", REAL_TYPE, PROBLEM, CELL_MODEL);
    initialize2DVariableFromFile(X_r1, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "X_r1", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastframeX_r2.txt", REAL_TYPE, PROBLEM, CELL_MODEL);
    initialize2DVariableFromFile(X_r2, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "X_r2", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastframeX_s.txt", REAL_TYPE, PROBLEM, CELL_MODEL);
    initialize2DVariableFromFile(X_s, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "X_s", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastframem.txt", REAL_TYPE, PROBLEM, CELL_MODEL);
    initialize2DVariableFromFile(m, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "m", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastframeh.txt", REAL_TYPE, PROBLEM, CELL_MODEL);
    initialize2DVariableFromFile(h, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "h", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastframej.txt", REAL_TYPE, PROBLEM, CELL_MODEL);
    initialize2DVariableFromFile(_j, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "j", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastframed.txt", REAL_TYPE, PROBLEM, CELL_MODEL);
    initialize2DVariableFromFile(d, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "d", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastframef.txt", REAL_TYPE, PROBLEM, CELL_MODEL);
    initialize2DVariableFromFile(f, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "f", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastframef2.txt", REAL_TYPE, PROBLEM, CELL_MODEL);
    initialize2DVariableFromFile(f2, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "f2", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastframefCaSS.txt", REAL_TYPE, PROBLEM, CELL_MODEL);
    initialize2DVariableFromFile(fCaSS, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "fCaSS", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastframes.txt", REAL_TYPE, PROBLEM, CELL_MODEL);
    initialize2DVariableFromFile(s, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "s", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastframer.txt", REAL_TYPE, PROBLEM, CELL_MODEL);
    initialize2DVariableFromFile(r, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "r", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastframeCa_i.txt", REAL_TYPE, PROBLEM, CELL_MODEL);
    initialize2DVariableFromFile(Ca_i, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "Ca_i", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastframeCa_SR.txt", REAL_TYPE, PROBLEM, CELL_MODEL);
    initialize2DVariableFromFile(Ca_SR, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "Ca_SR", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastframeCa_SS.txt", REAL_TYPE, PROBLEM, CELL_MODEL);
    initialize2DVariableFromFile(Ca_SS, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "Ca_SS", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastframeR_prime.txt", REAL_TYPE, PROBLEM, CELL_MODEL);
    initialize2DVariableFromFile(R_prime, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "R_prime", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastframeNa_i.txt", REAL_TYPE, PROBLEM, CELL_MODEL);
    initialize2DVariableFromFile(Na_i, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "Na_i", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastframeK_i.txt", REAL_TYPE, PROBLEM, CELL_MODEL);
    initialize2DVariableFromFile(K_i, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "K_i", real_ref_dx, real_def_dy);

#endif // TT2
    
        free(pathToRestoreStateFiles);

#endif // RESTORE_STATE

#ifdef SHIFT_STATE

    printf("\n");
    printf("Shifting state variables...\n");

    // Shift variables
    real lengthToShift = 0.5f;
    
    shift2DVariableToLeft(V, Nx, Ny, lengthToShift, delta_x, delta_y, V_init, "V");

#ifdef AFHN

    shift2DVariableToLeft(W, Nx, Ny, lengthToShift, delta_x, delta_y, W_init, "W");

#endif // AFHN
#ifdef TT2

    shift2DVariableToLeft(X_r1, Nx, Ny, lengthToShift, delta_x, delta_y, X_r1_init, "X_r1");
    shift2DVariableToLeft(X_r2, Nx, Ny, lengthToShift, delta_x, delta_y, X_r2_init, "X_r2");
    shift2DVariableToLeft(X_s, Nx, Ny, lengthToShift, delta_x, delta_y, X_s_init, "X_s");
    shift2DVariableToLeft(m, Nx, Ny, lengthToShift, delta_x, delta_y, m_init, "m");
    shift2DVariableToLeft(h, Nx, Ny, lengthToShift, delta_x, delta_y, h_init, "h");
    shift2DVariableToLeft(_j, Nx, Ny, lengthToShift, delta_x, delta_y, j_init, "j");
    shift2DVariableToLeft(d, Nx, Ny, lengthToShift, delta_x, delta_y, d_init, "d");
    shift2DVariableToLeft(f, Nx, Ny, lengthToShift, delta_x, delta_y, f_init, "f");
    shift2DVariableToLeft(f2, Nx, Ny, lengthToShift, delta_x, delta_y, f2_init, "f2");
    shift2DVariableToLeft(fCaSS, Nx, Ny, lengthToShift, delta_x, delta_y, fCaSS_init, "fCaSS");
    shift2DVariableToLeft(s, Nx, Ny, lengthToShift, delta_x, delta_y, s_init, "s");
    shift2DVariableToLeft(r, Nx, Ny, lengthToShift, delta_x, delta_y, r_init, "r");
    shift2DVariableToLeft(Ca_i, Nx, Ny, lengthToShift, delta_x, delta_y, Ca_i_init, "Ca_i");
    shift2DVariableToLeft(Ca_SR, Nx, Ny, lengthToShift, delta_x, delta_y, Ca_SR_init, "Ca_SR");
    shift2DVariableToLeft(Ca_SS, Nx, Ny, lengthToShift, delta_x, delta_y, Ca_SS_init, "Ca_SS");
    shift2DVariableToLeft(R_prime, Nx, Ny, lengthToShift, delta_x, delta_y, R_prime_init, "R_prime");
    shift2DVariableToLeft(Na_i, Nx, Ny, lengthToShift, delta_x, delta_y, Na_i_init, "Na_i");
    shift2DVariableToLeft(K_i, Nx, Ny, lengthToShift, delta_x, delta_y, K_i_init, "K_i");

#endif // TT2
#endif // SHIFT_STATE

    // Populate auxiliary arrays for Thomas algorithm
    real phi_x = delta_t / (delta_x * delta_x);
    real phi_y = delta_t / (delta_y * delta_y);
    real tau = 0.5f;

#ifndef TT2
    
    real diff_coeff = sigma / (Cm * chi);

#else // if def TT2

    real diff_coeff = sigma / chi;

#endif // not TT2

#if defined(SSIADI) || defined(THETASSIADI) || defined(OSADI)

    // Auxiliary arrays for Thomas algorithm
    real *la_x = (real *)malloc(Nx * sizeof(real)); // subdiagonal
    real *lb_x = (real *)malloc(Nx * sizeof(real)); // diagonal
    real *lc_x = (real *)malloc(Nx * sizeof(real)); // superdiagonal

    real *la_y = (real *)malloc(Ny * sizeof(real)); // subdiagonal
    real *lb_y = (real *)malloc(Ny * sizeof(real)); // diagonal
    real *lc_y = (real *)malloc(Ny * sizeof(real)); // superdiagonal

#if defined(SSIADI)

    populateDiagonalThomasAlgorithm(la_x, lb_x, lc_x, Nx, 0.5f * phi_x * diff_coeff);
    populateDiagonalThomasAlgorithm(la_y, lb_y, lc_y, Ny, 0.5f * phi_y * diff_coeff);

#endif // SSIADI

#if defined(THETASSIADI)

    populateDiagonalThomasAlgorithm(la_x, lb_x, lc_x, Nx, THETA * phi_x * diff_coeff);
    populateDiagonalThomasAlgorithm(la_y, lb_y, lc_y, Ny, THETA * phi_y * diff_coeff);
    tau = 1.0f - THETA;

#endif // THETASSIADI

#ifdef OSADI

    populateDiagonalThomasAlgorithm(la_x, lb_x, lc_x, Nx, phi_x * diff_coeff);
    populateDiagonalThomasAlgorithm(la_y, lb_y, lc_y, Ny, phi_y * diff_coeff);

#endif // OSADI

    // Prefactorization
    thomasFactorConstantBatch(la_x, lb_x, lc_x, Nx);
    thomasFactorConstantBatch(la_y, lb_y, lc_y, Ny);

    // Malloc and copy to device
    real *d_la_x, *d_lb_x, *d_lc_x, *d_la_y, *d_lb_y, *d_lc_y;
    CUDA_CALL(cudaMalloc(&d_la_x, Nx * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_lb_x, Nx * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_lc_x, Nx * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_la_y, Ny * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_lb_y, Ny * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_lc_y, Ny * sizeof(real)));

    CUDA_CALL(cudaMemcpy(d_la_x, la_x, Nx * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_lb_x, lb_x, Nx * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_lc_x, lc_x, Nx * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_la_y, la_y, Ny * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_lb_y, lb_y, Ny * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_lc_y, lc_y, Ny * sizeof(real), cudaMemcpyHostToDevice));

#endif // SSIADI || THETASSIADI || OSADI

    // Create device variables
    real *d_V, *d_RHS, *d_Vtilde, *d_partRHS;
    
    // Allocate memory on device
    CUDA_CALL(cudaMalloc(&d_V, Nx * Ny * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_RHS, Nx * Ny * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_Vtilde, Nx * Ny * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_partRHS, Nx * Ny * sizeof(real)));
    
    // Copy memory from host to device
    CUDA_CALL(cudaMemcpy(d_V, V, Nx * Ny * sizeof(real), cudaMemcpyHostToDevice));

#ifdef AFHN

    real *d_W;
    CUDA_CALL(cudaMalloc(&d_W, Nx * Ny * sizeof(real)));
    CUDA_CALL(cudaMemcpy(d_W, W, Nx * Ny * sizeof(real), cudaMemcpyHostToDevice));

#endif // AFHN

#ifdef TT2

    real *d_X_r1, *d_X_r2, *d_X_s, *d_m, *d_h, *d_j, *d_d, *d_f, *d_f2, *d_fCaSS, *d_s, *d_r, *d_Ca_i, *d_Ca_SR, *d_Ca_SS, *d_R_prime, *d_Na_i, *d_K_i;
    CUDA_CALL(cudaMalloc(&d_X_r1, Nx * Ny * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_X_r2, Nx * Ny * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_X_s, Nx * Ny * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_m, Nx * Ny * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_h, Nx * Ny * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_j, Nx * Ny * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_d, Nx * Ny * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_f, Nx * Ny * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_f2, Nx * Ny * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_fCaSS, Nx * Ny * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_s, Nx * Ny * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_r, Nx * Ny * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_Ca_i, Nx * Ny * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_Ca_SR, Nx * Ny * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_Ca_SS, Nx * Ny * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_R_prime, Nx * Ny * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_Na_i, Nx * Ny * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_K_i, Nx * Ny * sizeof(real)));
    CUDA_CALL(cudaMemcpy(d_X_r1, X_r1, Nx * Ny * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_X_r2, X_r2, Nx * Ny * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_X_s, X_s, Nx * Ny * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_m, m, Nx * Ny * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_h, h, Nx * Ny * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_j, _j, Nx * Ny * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_d, d, Nx * Ny * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_f, f, Nx * Ny * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_f2, f2, Nx * Ny * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_fCaSS, fCaSS, Nx * Ny * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_s, s, Nx * Ny * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_r, r, Nx * Ny * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_Ca_i, Ca_i, Nx * Ny * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_Ca_SR, Ca_SR, Nx * Ny * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_Ca_SS, Ca_SS, Nx * Ny * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_R_prime, R_prime, Nx * Ny * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_Na_i, Na_i, Nx * Ny * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_K_i, K_i, Nx * Ny * sizeof(real), cudaMemcpyHostToDevice));

#endif // TT2

#ifndef CONVERGENCE_ANALYSIS_FORCING_TERM
#ifdef MONODOMAIN

    // Allocate array for the stimuli
    Stimulus *stimuli = (Stimulus *)malloc(numberOfStimuli * sizeof(Stimulus));
    populateStimuli(stimuli, delta_x, delta_y);

    // Allocate and copy array for the stimuli on device
    Stimulus *d_stimuli;
    CUDA_CALL(cudaMalloc(&d_stimuli, numberOfStimuli * sizeof(Stimulus)));
    CUDA_CALL(cudaMemcpy(d_stimuli, stimuli, numberOfStimuli * sizeof(Stimulus), cudaMemcpyHostToDevice));

#endif // MONODOMAIN
#endif // not CONVERGENCE_ANALYSIS_FORCING_TERM

    // CUDA grid and block allocation
    // Device properties
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);

    // Number of SMs and minimum number of blocks to maximize the parallelism
    int numSMs = prop.multiProcessorCount;
    int minBlocks = 4 * numSMs;

    // Print information
    printf("\n");
    printf("Device name: %s (%d SMs)\n", prop.name, numSMs);

    // Calculate the number of blocks and threads for the full domain kernels
    dim3 blockSize(ceil(sqrt(BLOCK_SIZE)), ceil(sqrt(BLOCK_SIZE)));
    dim3 gridSize((Nx + blockSize.x - 1) / blockSize.x, (Ny + blockSize.y - 1) / blockSize.y);

    // Adjust the number of blocks
    if (gridSize.x * gridSize.y < minBlocks)
        gridSize.x = (minBlocks + gridSize.y - 1) / gridSize.y;
    
    // Print information
    printf("\n");
    printf("For full domain kernels:\n");
    printf("Block size: %d x %d threads (total %d threads per block)\n", blockSize.x, blockSize.y, blockSize.x * blockSize.y);
    printf("Grid size: %d x %d blocks (total %d blocks, total %d threads)\n", gridSize.x, gridSize.y, gridSize.x * gridSize.y, gridSize.x * gridSize.y * blockSize.x * blockSize.y);

#if defined(SSIADI) || defined(THETASSIADI) || defined(OSADI)

    // Calculate the number of blocks and threads for individual directions of ADI that will be used in Thomas kernel
    int gridSize_x = (Nx + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int gridSize_y = (Ny + BLOCK_SIZE - 1) / BLOCK_SIZE;

    // Adjust the number of blocks
    if (gridSize_x < minBlocks)
        gridSize_x = minBlocks;
    if (gridSize_y < minBlocks)
        gridSize_y = minBlocks;

    // Print information
    printf("\n");
    printf("For Thomas kernel:\n");
    printf("Grid size for x: %d blocks (%d threads per block, total %d threads)\n", gridSize_x, BLOCK_SIZE, gridSize_x * BLOCK_SIZE);
    printf("Grid size for y: %d blocks (%d threads per block, total %d threads)\n", gridSize_y, BLOCK_SIZE, gridSize_y * BLOCK_SIZE);

#endif // SSIADI || THETASSIADI || OSADI
    
    // Create directories
    char *pathToSaveData = (char *)malloc(MAX_STRING_SIZE * sizeof(char));
    createDirectories(delta_t, delta_x, delta_y, pathToSaveData);

#ifdef SAVE_FRAMES

    // Save frames variables
    char framesPath[MAX_STRING_SIZE];
    FILE *fpFrames;
    snprintf(framesPath, MAX_STRING_SIZE * sizeof(char), "%s/frames.txt", pathToSaveData);
    fpFrames = fopen(framesPath, "w");
    real startSaveFramesTime, finishSaveFramesTime, elapsedSaveFramesTime = 0.0f;

#endif // SAVE_FRAMES

#ifdef MEASURE_VELOCITY

    // Measure velocity
    real stim_velocity = 0.0;
    real first_point_time = 0.0;
    real last_point_time = 0.0;
    bool aux_stim_velocity_flag = false;
    bool stim_velocity_measured = false;
    real startMeasureVelocityTime, finishMeasureVelocityTime, elapsedMeasureVelocityTime = 0.0f;

#endif // MEASURE_VELOCITY

    // Variables for measuring the execution time
    real startExecutionTime, finishExecutionTime, elapsedExecutionTime = 0.0f;

    // Auxiliary variables for the time loop
    int timeStepCounter = 0;
    real actualTime = 0.0f;

    // Start measuring the execution time
    printf("\n");
    printf("Starting simulation...\n");
    startExecutionTime = omp_get_wtime();

#if defined(LINMONO) || defined(MONODOMAIN)

    while (timeStepCounter < M)
    {
        // Get time step
        actualTime = time[timeStepCounter];

        // ================================================!
        //  Calculate Approx. and ODEs                     !
        // ================================================!
        
#if defined(AFHN)
#ifdef CONVERGENCE_ANALYSIS_FORCING_TERM

        computeApprox<<<gridSize, blockSize>>>(Nx, Ny, delta_t, phi_x, phi_y, diff_coeff, delta_x, delta_y, actualTime, d_V, d_Vtilde, d_partRHS, d_W);

#else // not CONVERGENCE_ANALYSIS_FORCING_TERM

        computeApprox<<<gridSize, blockSize>>>(Nx, Ny, delta_t, phi_x, phi_y, diff_coeff, actualTime, d_V, d_partRHS, d_W, d_stimuli);

#endif // CONVERGENCE_ANALYSIS_FORCING_TERM
#endif // AFHN

#ifdef TT2

        computeApprox<<<gridSize, blockSize>>>(Nx, Ny, delta_t, phi_x, phi_y, diff_coeff, actualTime, d_V, d_partRHS, d_X_r1, d_X_r2, d_X_s, d_m, d_h, d_j, d_d, d_f, d_f2, d_fCaSS, d_s, d_r, d_Ca_i, d_Ca_SR, d_Ca_SS, d_R_prime, d_Na_i, d_K_i, d_stimuli);

#endif // TT2

        CUDA_CALL(cudaDeviceSynchronize());

#if defined(SSIADI) || defined(THETASSIADI)

        // ================================================!
        //  Calculate V at n+1/2 -> Result goes to RHS     !
        //  diffusion implicit in y and explicit in x      !
        // ================================================!
        prepareRHSwithjDiff<<<gridSize, blockSize>>>(Nx, Ny, phi_x, diff_coeff, tau, d_V, d_RHS, d_partRHS);
        CUDA_CALL(cudaDeviceSynchronize());

        parallelThomas<<<gridSize_y, BLOCK_SIZE>>>(Nx, Ny, d_RHS, d_la_y, d_lb_y, d_lc_y);
        CUDA_CALL(cudaDeviceSynchronize());

#endif // SSIADI || THETASSIADI

        // =========================================================!
        //  Transpose d_RHS to d_Vtilde (d_Vtilde as an aux var)    !
        // =========================================================!
        parallelTranspose<<<gridSize, blockSize>>>(Nx, Ny, d_RHS, d_Vtilde);
        CUDA_CALL(cudaDeviceSynchronize());

#if defined(SSIADI) || defined(THETASSIADI)

        // ================================================!
        //  Calculate V at n+1 -> Result goes to V         !
        //  diffusion implicit in x and explicit in y      !
        // ================================================!
        prepareRHSwithiDiff<<<gridSize, blockSize>>>(Nx, Ny, phi_y, diff_coeff, tau, d_Vtilde, d_V, d_partRHS);
        CUDA_CALL(cudaDeviceSynchronize());

        parallelThomas<<<gridSize_x, BLOCK_SIZE>>>(Ny, Nx, d_V, d_la_x, d_lb_x, d_lc_x);
        CUDA_CALL(cudaDeviceSynchronize());
        
#endif // SSIADI || THETASSIADI

#ifdef SAVE_FRAMES

        startSaveFramesTime = omp_get_wtime();

        // If save frames is true and time step is multiple of frame save rate
        if (timeStepCounter % frameSaveRate == 0)
        {
            // Copy memory of d_V from device to host V
            CUDA_CALL(cudaMemcpy(V, d_V, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));
            CUDA_CALL(cudaDeviceSynchronize());

            // Save frame
            fprintf(fpFrames, "%lf\n", actualTime);
            saveFrame(fpFrames, V, Nx, Ny);
            SUCCESSMSG("Frame at time %.2f ms saved to %s\n", actualTime, framesPath);
        }

        finishSaveFramesTime = omp_get_wtime();
        elapsedSaveFramesTime += finishSaveFramesTime - startSaveFramesTime;
        
#endif // SAVE_FRAMES

#ifdef MEASURE_VELOCITY

        startMeasureVelocityTime = omp_get_wtime();

        // Calculate stim velocity
        if (!stim_velocity_measured)
        {
            // Copy memory of d_V from device to host V
            CUDA_CALL(cudaMemcpy(V, d_V, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));
            CUDA_CALL(cudaDeviceSynchronize());

            real begin = Lx / 3.0f;
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
                    INFOMSG("Stim velocity (measured from %.2f to %.2f cm) is %.5g m/s\n", begin, end, stim_velocity);
                }
            }
        }

        finishMeasureVelocityTime = omp_get_wtime();
        elapsedMeasureVelocityTime += finishMeasureVelocityTime - startMeasureVelocityTime;

#endif // MEASURE_VELOCITY

        // Update time step counter
        timeStepCounter++;
    }
    
#endif // LINMONO || MONODOMAIN

    finishExecutionTime = omp_get_wtime();
    elapsedExecutionTime += finishExecutionTime - startExecutionTime;

    // Copy memory of d_V from device to host V
    CUDA_CALL(cudaMemcpy(V, d_V, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));

#ifdef SAVE_FRAMES

    fprintf(fpFrames, "%lf\n", actualTime);
    saveFrame(fpFrames, V, Nx, Ny);
    SUCCESSMSG("Frame at time %.2f ms saved to %s\n", actualTime, framesPath);
    fclose(fpFrames);

#endif // SAVE_FRAMES

    printf("Simulation done!\n");
    printf("\n");

#ifdef CONVERGENCE_ANALYSIS_FORCING_TERM

    real **exact = (real **)malloc(Ny * sizeof(real *));
    for (int i = 0; i < Ny; i++)
        exact[i] = (real *)malloc(Nx * sizeof(real));

    real norm2error = calculateNorm2Error(V, exact, Nx, Ny, totalTime, delta_x, delta_y);

    char exactFilePath[MAX_STRING_SIZE];
    snprintf(exactFilePath, MAX_STRING_SIZE * sizeof(char), "%s/exact.txt", pathToSaveData);
    FILE *fpExact = fopen(exactFilePath, "w");
    
    char errorsFilePath[MAX_STRING_SIZE];
    snprintf(errorsFilePath, MAX_STRING_SIZE * sizeof(char), "%s/errors.txt", pathToSaveData);
    FILE *fpErrors = fopen(errorsFilePath, "w");

    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            int index = i * Nx + j;
            fprintf(fpExact, "%e ", exact[i][j]);
            fprintf(fpErrors, "%e ", abs(V[index] - exact[i][j]));
        }
        fprintf(fpExact, "\n");
        fprintf(fpErrors, "\n");
    }

    for (int i = 0; i < Ny; i++)
        free(exact[i]);
    free(exact);

    SUCCESSMSG("Exact solution saved to %s\n", exactFilePath);
    SUCCESSMSG("Errors saved to %s\n", errorsFilePath);
    fclose(fpExact);
    fclose(fpErrors);

#endif // CONVERGENCE_ANALYSIS_FORCING_TERM

    // Write infos to file
    char infosFilePath[MAX_STRING_SIZE];
    snprintf(infosFilePath, MAX_STRING_SIZE * sizeof(char), "%s/infos.txt", pathToSaveData);
    FILE *fpInfos = fopen(infosFilePath, "w");
    fprintf(fpInfos, "EXECUTION TYPE = %s\n", EXECUTION_TYPE);
    fprintf(fpInfos, "PRECISION = %s\n", REAL_TYPE);
    fprintf(fpInfos, "PROBLEM = %s\n", PROBLEM);
    fprintf(fpInfos, "CELL MODEL = %s\n", CELL_MODEL);
    fprintf(fpInfos, "METHOD = %s\n", METHOD);

#ifdef THETA

    fprintf(fpInfos, "theta = %.2f\n", THETA);
    
#endif // THETA

    fprintf(fpInfos, "\n");
    fprintf(fpInfos, "DOMAIN LENGTH IN X = %.4g cm\n", Lx);
    fprintf(fpInfos, "DOMAIN LENGTH IN Y = %.4g cm\n", Ly);
    fprintf(fpInfos, "TOTAL TIME = %.4g ms\n", totalTime);

    fprintf(fpInfos, "\n");
    fprintf(fpInfos, "delta_t = %.5g ms (%d time steps)\n", delta_t, M);
    fprintf(fpInfos, "delta_x = %.5g cm (%d um) (%d space steps in x)\n", delta_x, CM_TO_UM(delta_x), Nx);
    fprintf(fpInfos, "delta_y = %.5g cm (%d um) (%d space steps in y)\n", delta_y, CM_TO_UM(delta_y), Ny);
    fprintf(fpInfos, "TOTAL POINTS IN DOMAIN = %d\n", Nx * Ny);  

    fprintf(fpInfos, "\n");
    fprintf(fpInfos, "DEVICE NAME = %s (%d SMs)\n", prop.name, numSMs);

    fprintf(fpInfos, "\n");
    fprintf(fpInfos, "FOR FULL DOMAIN KERNELS:\n");
    fprintf(fpInfos, "BLOCK SIZE = %d x %d threads (%d threads per block)\n", blockSize.x, blockSize.y, blockSize.x * blockSize.y);
    fprintf(fpInfos, "GRID SIZE = %d x %d blocks (total %d blocks, total %d threads)\n", gridSize.x, gridSize.y, gridSize.x * gridSize.y, gridSize.x * gridSize.y * blockSize.x * blockSize.y);

#if defined(SSIADI) || defined(THETASSIADI) || defined(OSADI)

    fprintf(fpInfos, "\n");
    fprintf(fpInfos, "FOR THOMAS KERNEL:\n");
    fprintf(fpInfos, "GRID SIZE FOR X = %d blocks (%d threads per block, total %d threads)\n", gridSize_x, BLOCK_SIZE, gridSize_x * BLOCK_SIZE);
    fprintf(fpInfos, "GRID SIZE FOR Y = %d blocks (%d threads per block, total %d threads)\n", gridSize_y, BLOCK_SIZE, gridSize_y * BLOCK_SIZE);

#endif // SSIADI || THETASSIADI || OSADI

#ifdef MEASURE_VELOCITY

    fprintf(fpInfos, "\n");
    fprintf(fpInfos, "STIMULUS VELOCITY = %.5g m/s\n", stim_velocity);

#endif // MEASURE_VELOCITY

#ifdef CONVERGENCE_ANALYSIS_FORCING_TERM

    fprintf(fpInfos, "\n");
    fprintf(fpInfos, "NORM-2 ERROR = %lf\n", norm2error);

#endif // CONVERGENCE_ANALYSIS_FORCING_TERM
    
    fprintf(fpInfos, "\n");

#ifdef MEASURE_VELOCITY

    fprintf(fpInfos, "TIME TO MEASURE VELOCITY = %.5g s\n", elapsedMeasureVelocityTime);

#endif // MEASURE_VELOCITY

#ifdef SAVE_FRAMES

    fprintf(fpInfos, "TIME TO SAVE FRAMES = %.5g s\n", elapsedSaveFramesTime);

#endif // SAVE_FRAMES

    fprintf(fpInfos, "SIMULATION TOTAL EXECUTION TIME = %.5g s\n", elapsedExecutionTime);
    fclose(fpInfos);

    INFOMSG("Simulation total execution time = %.5g s\n", elapsedExecutionTime);
    SUCCESSMSG("Simulation infos saved to %s\n", infosFilePath);

#ifdef SAVE_LAST_FRAME

    // Save last frame
    char lastFrameFilePath[MAX_STRING_SIZE];
    snprintf(lastFrameFilePath, MAX_STRING_SIZE * sizeof(char), "%s/lastframe.txt", pathToSaveData);
    FILE *fpLast = fopen(lastFrameFilePath, "w");

    saveFrame(fpLast, V, Nx, Ny);

    SUCCESSMSG("Last frame saved to %s\n", lastFrameFilePath);
    fclose(fpLast);

#endif // SAVE_LAST_FRAME

#ifdef SAVE_LAST_STATE
#ifdef AFHN

    char lastFrameFilePathV[MAX_STRING_SIZE], lastFrameFilePathW[MAX_STRING_SIZE];
    snprintf(lastFrameFilePathV, MAX_STRING_SIZE * sizeof(char), "%s/lastframeV.txt", pathToSaveData);
    snprintf(lastFrameFilePathW, MAX_STRING_SIZE * sizeof(char), "%s/lastframeW.txt", pathToSaveData);

    CUDA_CALL(cudaMemcpy(W, d_W, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaFree(d_W));

    FILE *fpLastV = fopen(lastFrameFilePathV, "w");
    FILE *fpLastW = fopen(lastFrameFilePathW, "w");

    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            int index = i * Nx + j;
            fprintf(fpLastV, "%e ", V[index]);
            fprintf(fpLastW, "%e ", W[index]);
        }
        fprintf(fpLastV, "\n");
        fprintf(fpLastW, "\n");
    }

    SUCCESSMSG("Last V frame saved to %s\n", lastFrameFilePathV);
    SUCCESSMSG("Last W frame saved to %s\n", lastFrameFilePathW);
    fclose(fpLastV);
    fclose(fpLastW);

#endif // AFHN

#ifdef TT2

    char lastFrameFilePathV[MAX_STRING_SIZE], lastFrameFilePathX_r1[MAX_STRING_SIZE], lastFrameFilePathX_r2[MAX_STRING_SIZE], lastFrameFilePathX_s[MAX_STRING_SIZE],
         lastFrameFilePathm[MAX_STRING_SIZE], lastFrameFilePathh[MAX_STRING_SIZE], lastFrameFilePathj[MAX_STRING_SIZE], lastFrameFilePathd[MAX_STRING_SIZE], lastFrameFilePathf[MAX_STRING_SIZE],
         lastFrameFilePathf2[MAX_STRING_SIZE], lastFrameFilePathfCaSS[MAX_STRING_SIZE], lastFrameFilePaths[MAX_STRING_SIZE], lastFrameFilePathr[MAX_STRING_SIZE],
         lastFrameFilePathR_prime[MAX_STRING_SIZE], lastFrameFilePathCa_i[MAX_STRING_SIZE], lastFrameFilePathCa_SR[MAX_STRING_SIZE], lastFrameFilePathCa_SS[MAX_STRING_SIZE],
         lastFrameFilePathNa_i[MAX_STRING_SIZE], lastFrameFilePathK_i[MAX_STRING_SIZE];
    snprintf(lastFrameFilePathV, MAX_STRING_SIZE * sizeof(char), "%s/lastframeV.txt", pathToSaveData);
    snprintf(lastFrameFilePathX_r1, MAX_STRING_SIZE * sizeof(char), "%s/lastframeX_r1.txt", pathToSaveData);
    snprintf(lastFrameFilePathX_r2, MAX_STRING_SIZE * sizeof(char), "%s/lastframeX_r2.txt", pathToSaveData);
    snprintf(lastFrameFilePathX_s, MAX_STRING_SIZE * sizeof(char), "%s/lastframeX_s.txt", pathToSaveData);
    snprintf(lastFrameFilePathm, MAX_STRING_SIZE * sizeof(char), "%s/lastframem.txt", pathToSaveData);
    snprintf(lastFrameFilePathh, MAX_STRING_SIZE * sizeof(char), "%s/lastframeh.txt", pathToSaveData);
    snprintf(lastFrameFilePathj, MAX_STRING_SIZE * sizeof(char), "%s/lastframej.txt", pathToSaveData);
    snprintf(lastFrameFilePathd, MAX_STRING_SIZE * sizeof(char), "%s/lastframed.txt", pathToSaveData);
    snprintf(lastFrameFilePathf, MAX_STRING_SIZE * sizeof(char), "%s/lastframef.txt", pathToSaveData);
    snprintf(lastFrameFilePathf2, MAX_STRING_SIZE * sizeof(char), "%s/lastframef2.txt", pathToSaveData);
    snprintf(lastFrameFilePathfCaSS, MAX_STRING_SIZE * sizeof(char), "%s/lastframefCaSS.txt", pathToSaveData);
    snprintf(lastFrameFilePaths, MAX_STRING_SIZE * sizeof(char), "%s/lastframes.txt", pathToSaveData);
    snprintf(lastFrameFilePathr, MAX_STRING_SIZE * sizeof(char), "%s/lastframer.txt", pathToSaveData);
    snprintf(lastFrameFilePathR_prime, MAX_STRING_SIZE * sizeof(char), "%s/lastframeR_prime.txt", pathToSaveData);
    snprintf(lastFrameFilePathCa_i, MAX_STRING_SIZE * sizeof(char), "%s/lastframeCa_i.txt", pathToSaveData);
    snprintf(lastFrameFilePathCa_SR, MAX_STRING_SIZE * sizeof(char), "%s/lastframeCa_SR.txt", pathToSaveData);
    snprintf(lastFrameFilePathCa_SS, MAX_STRING_SIZE * sizeof(char), "%s/lastframeCa_SS.txt", pathToSaveData);
    snprintf(lastFrameFilePathNa_i, MAX_STRING_SIZE * sizeof(char), "%s/lastframeNa_i.txt", pathToSaveData);
    snprintf(lastFrameFilePathK_i, MAX_STRING_SIZE * sizeof(char), "%s/lastframeK_i.txt", pathToSaveData);

    CUDA_CALL(cudaMemcpy(X_r1, d_X_r1, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(X_r2, d_X_r2, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(X_s, d_X_s, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(m, d_m, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(h, d_h, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(_j, d_j, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(d, d_d, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(f, d_f, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(f2, d_f2, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(fCaSS, d_fCaSS, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(s, d_s, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(r, d_r, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(R_prime, d_R_prime, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(Ca_i, d_Ca_i, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(Ca_SR, d_Ca_SR, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(Ca_SS, d_Ca_SS, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(Na_i, d_Na_i, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(K_i, d_K_i, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));

    FILE *fpLastV = fopen(lastFrameFilePathV, "w");
    FILE *fpLastX_r1 = fopen(lastFrameFilePathX_r1, "w");
    FILE *fpLastX_r2 = fopen(lastFrameFilePathX_r2, "w");
    FILE *fpLastX_s = fopen(lastFrameFilePathX_s, "w");
    FILE *fpLastm = fopen(lastFrameFilePathm, "w");
    FILE *fpLasth = fopen(lastFrameFilePathh, "w");
    FILE *fpLastj = fopen(lastFrameFilePathj, "w");
    FILE *fpLastd = fopen(lastFrameFilePathd, "w");
    FILE *fpLastf = fopen(lastFrameFilePathf, "w");
    FILE *fpLastf2 = fopen(lastFrameFilePathf2, "w");
    FILE *fpLastfCaSS = fopen(lastFrameFilePathfCaSS, "w");
    FILE *fpLasts = fopen(lastFrameFilePaths, "w");
    FILE *fpLastr = fopen(lastFrameFilePathr, "w");
    FILE *fpLastR_prime = fopen(lastFrameFilePathR_prime, "w");
    FILE *fpLastCa_i = fopen(lastFrameFilePathCa_i, "w");
    FILE *fpLastCa_SR = fopen(lastFrameFilePathCa_SR, "w");
    FILE *fpLastCa_SS = fopen(lastFrameFilePathCa_SS, "w");
    FILE *fpLastNa_i = fopen(lastFrameFilePathNa_i, "w");
    FILE *fpLastK_i = fopen(lastFrameFilePathK_i, "w");

    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            int index = i * Nx + j;
            fprintf(fpLastV, "%e ", V[index]);
            fprintf(fpLastX_r1, "%e ", X_r1[index]);
            fprintf(fpLastX_r2, "%e ", X_r2[index]);
            fprintf(fpLastX_s, "%e ", X_s[index]);
            fprintf(fpLastm, "%e ", m[index]);
            fprintf(fpLasth, "%e ", h[index]);
            fprintf(fpLastj, "%e ", _j[index]);
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
        }
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
    }

    SUCCESSMSG("Last V frame saved to %s\n", lastFrameFilePathV);
    SUCCESSMSG("Last X_r1 frame saved to %s\n", lastFrameFilePathX_r1);
    SUCCESSMSG("Last X_r2 frame saved to %s\n", lastFrameFilePathX_r2);
    SUCCESSMSG("Last X_s frame saved to %s\n", lastFrameFilePathX_s);
    SUCCESSMSG("Last m frame saved to %s\n", lastFrameFilePathm);
    SUCCESSMSG("Last h frame saved to %s\n", lastFrameFilePathh);
    SUCCESSMSG("Last j frame saved to %s\n", lastFrameFilePathj);
    SUCCESSMSG("Last d frame saved to %s\n", lastFrameFilePathd);
    SUCCESSMSG("Last f frame saved to %s\n", lastFrameFilePathf);
    SUCCESSMSG("Last f2 frame saved to %s\n", lastFrameFilePathf2);
    SUCCESSMSG("Last fCaSS frame saved to %s\n", lastFrameFilePathfCaSS);
    SUCCESSMSG("Last s frame saved to %s\n", lastFrameFilePaths);
    SUCCESSMSG("Last r frame saved to %s\n", lastFrameFilePathr);
    SUCCESSMSG("Last R_prime frame saved to %s\n", lastFrameFilePathR_prime);
    SUCCESSMSG("Last Ca_i frame saved to %s\n", lastFrameFilePathCa_i);
    SUCCESSMSG("Last Ca_SR frame saved to %s\n", lastFrameFilePathCa_SR);
    SUCCESSMSG("Last Ca_SS frame saved to %s\n", lastFrameFilePathCa_SS);
    SUCCESSMSG("Last Na_i frame saved to %s\n", lastFrameFilePathNa_i);
    SUCCESSMSG("Last K_i frame saved to %s\n", lastFrameFilePathK_i);

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
    free(pathToSaveData);

    free(V);
    free(Vtilde);
    free(RHS);
    free(partRHS);
    CUDA_CALL(cudaFree(d_V));
    CUDA_CALL(cudaFree(d_Vtilde));
    CUDA_CALL(cudaFree(d_RHS));
    CUDA_CALL(cudaFree(d_partRHS));

#if defined(SSIADI) || defined(THETASSIADI) || defined(OSADI)

    free(la_x);
    free(lb_x);
    free(lc_x);
    free(la_y);
    free(lb_y);
    free(lc_y);
    CUDA_CALL(cudaFree(d_la_x));
    CUDA_CALL(cudaFree(d_lb_x));
    CUDA_CALL(cudaFree(d_lc_x));
    CUDA_CALL(cudaFree(d_la_y));
    CUDA_CALL(cudaFree(d_lb_y));
    CUDA_CALL(cudaFree(d_lc_y));

#endif // SSIADI || THETASSIADI || OSADI

#ifdef MONODOMAIN
#ifndef CONVERGENCE_ANALYSIS_FORCING_TERM

    free(stimuli);
    CUDA_CALL(cudaFree(d_stimuli));

#endif // not CONVERGENCE_ANALYSIS_FORCING_TERM

#ifdef AFHN

    free(W);
    CUDA_CALL(cudaFree(d_W));

#endif // AFHN

#ifdef TT2

    free(X_r1);
    free(X_r2);
    free(X_s);
    free(m);
    free(h);
    free(_j);
    free(d);
    free(f);
    free(f2);
    free(fCaSS);
    free(s);
    free(r);
    free(R_prime);
    free(Ca_i);
    free(Ca_SR);
    free(Ca_SS);
    free(Na_i);
    free(K_i);
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
#endif // MONODOMAIN

    // Reset device
    CUDA_CALL(cudaDeviceReset());

    return;
}

#endif // GPU_METHODS_H