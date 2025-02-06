#ifndef GPU_METHODS_H
#define GPU_METHODS_H

#include "auxfuncs.h"

void runSimulationGPU(real delta_t, real delta_x, real delta_y)
{
    // Number of steps
    int M = round(totalTime / delta_t); // Number of time steps
    int Nx = round(Lx / delta_x) + 1;   // Spatial steps in x
    int Ny = round(Ly / delta_y) + 1;   // Spatial steps in y

    printf("Nx = %d\n", Nx);
    printf("Ny = %d\n", Ny);
    printf("Points in the domain = %d\n", Nx * Ny);

    // Allocate and populate time array
    real *time = (real *)malloc(M * sizeof(real));
    initializeTimeArray(time, M, delta_t);

    // Allocate arrays for variables (2D matrices will be flattened) and populate them

#ifdef AFHN

    real *Vm = (real *)malloc(Nx * Ny * sizeof(real));
    real *W = (real *)malloc(Nx * Ny * sizeof(real));

    initialize2DVariableWithValue(Vm, Nx, Ny, Vm_init);
    initialize2DVariableWithValue(W, Nx, Ny, W_init);

#endif // AFHN

#ifdef TT2

    real *Vm = (real *)malloc(Nx * Ny * sizeof(real));
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

    initialize2DVariableWithValue(Vm, Nx, Ny, Vm_init);
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

#ifdef MV

    real *Vm = (real *)malloc(Nx * Ny * sizeof(real));
    real *v = (real *)malloc(Nx * Ny * sizeof(real));
    real *w = (real *)malloc(Nx * Ny * sizeof(real));
    real *s = (real *)malloc(Nx * Ny * sizeof(real));

    initialize2DVariableWithValue(Vm, Nx, Ny, u_init);
    initialize2DVariableWithValue(v, Nx, Ny, v_init);
    initialize2DVariableWithValue(w, Nx, Ny, w_init);
    initialize2DVariableWithValue(s, Nx, Ny, s_init);

#endif // MV

#ifdef RESTORE_STATE

    printf("\n");
    printf("Restoring state variables...\n");

    // Initialize variables with a solution
    real real_ref_dx = 0.0005f;
    real real_def_dy = 0.0005f;

    char *pathToRestoreStateFiles = (char *)malloc(MAX_STRING_SIZE * sizeof(char));

#ifdef AFHN

    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastframeVm.txt", REAL_TYPE, PROBLEM, CELL_MODEL);
    initialize2DVariableFromFile(Vm, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "Vm", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastframeW.txt", REAL_TYPE, PROBLEM, CELL_MODEL);
    initialize2DVariableFromFile(W, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "W", real_ref_dx, real_def_dy);

#endif // AFHN

#ifdef TT2

    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastframeVm.txt", REAL_TYPE, PROBLEM, CELL_MODEL);
    initialize2DVariableFromFile(Vm, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "Vm", real_ref_dx, real_def_dy);
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

#ifdef MV

    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastframeVm.txt", REAL_TYPE, PROBLEM, CELL_MODEL);
    initialize2DVariableFromFile(Vm, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "Vm", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastframev.txt", REAL_TYPE, PROBLEM, CELL_MODEL);
    initialize2DVariableFromFile(v, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "v", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastframew.txt", REAL_TYPE, PROBLEM, CELL_MODEL);
    initialize2DVariableFromFile(w, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "w", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastframes.txt", REAL_TYPE, PROBLEM, CELL_MODEL);
    initialize2DVariableFromFile(s, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "s", real_ref_dx, real_def_dy);

#endif // MV

    free(pathToRestoreStateFiles);

#endif // RESTORE_STATE

#ifdef SHIFT_STATE

    printf("\n");
    printf("Shifting state variables...\n");

    // Shift variables
    real lengthToShift = 0.5f;

#ifdef AFHN

    shift2DVariableToLeft(Vm, Nx, Ny, lengthToShift, delta_x, delta_y, Vm_init, "Vm");
    shift2DVariableToLeft(W, Nx, Ny, lengthToShift, delta_x, delta_y, W_init, "W");

#endif // AFHN

#ifdef TT2

    shift2DVariableToLeft(Vm, Nx, Ny, lengthToShift, delta_x, delta_y, Vm_init, "Vm");
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

#ifdef MV

    shift2DVariableToLeft(Vm, Nx, Ny, lengthToShift, delta_x, delta_y, u_init, "Vm");
    shift2DVariableToLeft(v, Nx, Ny, lengthToShift, delta_x, delta_y, v_init, "v");
    shift2DVariableToLeft(w, Nx, Ny, lengthToShift, delta_x, delta_y, w_init, "w");
    shift2DVariableToLeft(s, Nx, Ny, lengthToShift, delta_x, delta_y, s_init, "s");

#endif // MV

#endif // SHIFT_STATE

    // Populate auxiliary arrays for Thomas algorithm
    real phi_x = delta_t / (delta_x * delta_x);
    real phi_y = delta_t / (delta_y * delta_y);
    real tau = 0.5f; // Used for calculating the explicit diffusion term on the right-hand side of the ADI method

#ifdef AFHN

    real diff_coeff = sigma / (Cm * chi);

#endif // AFHN

#ifdef TT2

    real diff_coeff = sigma / chi;

#endif // TT2

#ifdef MV

    real diff_coeff = Dtilde / chi;

#endif // MV

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

    // Prefactorize the arrays for Thomas algorithm
    real *c_prime_x = (real *)malloc(Nx * sizeof(real));
    real *c_prime_y = (real *)malloc(Ny * sizeof(real));
    real *denominator_x = (real *)malloc(Nx * sizeof(real));
    real *denominator_y = (real *)malloc(Ny * sizeof(real));
    prefactorizeThomas(la_x, lb_x, lc_x, c_prime_x, denominator_x, Nx);
    prefactorizeThomas(la_y, lb_y, lc_y, c_prime_y, denominator_y, Ny);

    real *d_c_prime_x, *d_c_prime_y, *d_denominator_x, *d_denominator_y;
    CUDA_CALL(cudaMalloc(&d_c_prime_x, Nx * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_c_prime_y, Ny * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_denominator_x, Nx * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_denominator_y, Ny * sizeof(real)));

    CUDA_CALL(cudaMemcpy(d_c_prime_x, c_prime_x, Nx * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_c_prime_y, c_prime_y, Ny * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_denominator_x, denominator_x, Nx * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_denominator_y, denominator_y, Ny * sizeof(real), cudaMemcpyHostToDevice));

    free(c_prime_x);
    free(c_prime_y);
    free(denominator_x);
    free(denominator_y);

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

    // Free host memory
    free(la_x);
    free(lb_x);
    free(lc_x);
    free(la_y);
    free(lb_y);
    free(lc_y);

#endif // SSIADI || THETASSIADI || OSADI

    // Create device variable and allocate memory on device
    real *d_partRHS;
    CUDA_CALL(cudaMalloc(&d_partRHS, Nx * Ny * sizeof(real)));

#if defined(SSIADI) || defined(THETASSIADI)

    // Allocate memory for the RHS
    real *d_RHS;
    CUDA_CALL(cudaMalloc(&d_RHS, Nx * Ny * sizeof(real)));

#endif // SSIADI || THETASSIADI

#ifdef AFHN

    real *d_Vm, *d_W;
    CUDA_CALL(cudaMalloc(&d_Vm, Nx * Ny * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_W, Nx * Ny * sizeof(real)));

    // Copy data to device
    CUDA_CALL(cudaMemcpy(d_Vm, Vm, Nx * Ny * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_W, W, Nx * Ny * sizeof(real), cudaMemcpyHostToDevice));

#endif // AFHN

#ifdef TT2

    real *d_Vm, *d_X_r1, *d_X_r2, *d_X_s, *d_m, *d_h, *d_j, *d_d, *d_f, *d_f2, *d_fCaSS, *d_s, *d_r, *d_Ca_i, *d_Ca_SR, *d_Ca_SS, *d_R_prime, *d_Na_i, *d_K_i;
    CUDA_CALL(cudaMalloc(&d_Vm, Nx * Ny * sizeof(real)));
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

    // Copy data to device
    CUDA_CALL(cudaMemcpy(d_Vm, Vm, Nx * Ny * sizeof(real), cudaMemcpyHostToDevice));
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

#ifdef MV

    real *d_Vm, *d_v, *d_w, *d_s;
    CUDA_CALL(cudaMalloc(&d_Vm, Nx * Ny * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_v, Nx * Ny * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_w, Nx * Ny * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_s, Nx * Ny * sizeof(real)));

    // Copy data to device
    CUDA_CALL(cudaMemcpy(d_Vm, Vm, Nx * Ny * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_v, v, Nx * Ny * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_w, w, Nx * Ny * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_s, s, Nx * Ny * sizeof(real), cudaMemcpyHostToDevice));

#endif // MV

#ifdef MONODOMAIN

    // Allocate array for the stimuli
    Stimulus *stimuli = (Stimulus *)malloc(numberOfStimuli * sizeof(Stimulus));
    populateStimuli(stimuli, delta_x, delta_y);

    // Allocate and copy array for the stimuli on device
    Stimulus *d_stimuli;
    CUDA_CALL(cudaMalloc(&d_stimuli, numberOfStimuli * sizeof(Stimulus)));
    CUDA_CALL(cudaMemcpy(d_stimuli, stimuli, numberOfStimuli * sizeof(Stimulus), cudaMemcpyHostToDevice));

#endif // MONODOMAIN

    // CUDA grid and block allocation
    // Device properties
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);

    // Number of SMs and minimum number of blocks to maximize the parallelism
    int numSMs = prop.multiProcessorCount;
    int minBlocks = 2 * numSMs;

    // Print information
    printf("\n");
    printf("Device name: %s (%d SMs)\n", prop.name, numSMs);

    // Calculate the number of blocks and threads for the full domain kernels
    dim3 fullDomainBlockSize(FULL_DOMAIN_BLOCK_SIZE_X, FULL_DOMAIN_BLOCK_SIZE_Y);
    dim3 fullDomainGridSize((Nx + fullDomainBlockSize.x - 1) / fullDomainBlockSize.x, (Ny + fullDomainBlockSize.y - 1) / fullDomainBlockSize.y);

    // Adjust the number of blocks
    if (fullDomainGridSize.x * fullDomainGridSize.y < minBlocks)
        fullDomainGridSize.x = (minBlocks + fullDomainGridSize.y - 1) / fullDomainGridSize.y;

    // Print information
    printf("\n");
    printf("For full domain kernels:\n");
    printf("Block size: %d x %d threads (total %d threads per block)\n", fullDomainBlockSize.x, fullDomainBlockSize.y, fullDomainBlockSize.x * fullDomainBlockSize.y);
    printf("Grid size: %d x %d blocks (total %d blocks, total %d threads)\n", fullDomainGridSize.x, fullDomainGridSize.y, fullDomainGridSize.x * fullDomainGridSize.y, fullDomainGridSize.x * fullDomainGridSize.y * fullDomainBlockSize.x * fullDomainBlockSize.y);

#if defined(SSIADI) || defined(THETASSIADI) || defined(OSADI)

    // Calculate the number of blocks and threads for individual directions of ADI that will be used in Thomas kernel
    int gridSize_x = (Nx + THOMAS_KERNEL_BLOCK_SIZE - 1) / THOMAS_KERNEL_BLOCK_SIZE;
    int gridSize_y = (Ny + THOMAS_KERNEL_BLOCK_SIZE - 1) / THOMAS_KERNEL_BLOCK_SIZE;

    // Print information
    printf("\n");
    printf("For Thomas kernel:\n");
    printf("Grid size for x: %d blocks (%d threads per block, total %d threads)\n", gridSize_x, THOMAS_KERNEL_BLOCK_SIZE, gridSize_x * THOMAS_KERNEL_BLOCK_SIZE);
    printf("Grid size for y: %d blocks (%d threads per block, total %d threads)\n", gridSize_y, THOMAS_KERNEL_BLOCK_SIZE, gridSize_y * THOMAS_KERNEL_BLOCK_SIZE);

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
    real begin_point = Lx / 3.0f;
    real end_point = 2.0f * begin_point;
    int begin_point_index = round(begin_point / delta_x) + 1;
    int end_point_index = round(end_point / delta_x) + 1;
    real begin_point_time = 0.0;
    real end_point_time = 0.0;
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

#if defined(MONODOMAIN)

    while (timeStepCounter < M)
    {
        // Get time step
        actualTime = time[timeStepCounter];

        // ================================================!
        //  Calculate Approx. and ODEs                     !
        // ================================================!

#ifdef AFHN

        computeApprox<<<fullDomainGridSize, fullDomainBlockSize>>>(Nx, Ny, delta_t, phi_x, phi_y, diff_coeff, actualTime, d_stimuli, d_Vm, d_partRHS, d_W);

#endif // AFHN

#ifdef TT2

        computeApprox<<<fullDomainGridSize, fullDomainBlockSize>>>(Nx, Ny, delta_t, phi_x, phi_y, diff_coeff, actualTime, d_stimuli, d_Vm, d_partRHS, d_X_r1, d_X_r2, d_X_s, d_m, d_h, d_j, d_d, d_f, d_f2, d_fCaSS, d_s, d_r, d_Ca_i, d_Ca_SR, d_Ca_SS, d_R_prime, d_Na_i, d_K_i);

#endif // TT2

#ifdef MV

        computeApprox<<<fullDomainGridSize, fullDomainBlockSize>>>(Nx, Ny, delta_t, phi_x, phi_y, diff_coeff, actualTime, d_stimuli, d_Vm, d_partRHS, d_v, d_w, d_s);

#endif // MV

        CUDA_CALL(cudaDeviceSynchronize());

#if defined(SSIADI) || defined(THETASSIADI)

        // ================================================!
        //  Calculate Vm at n+1/2 -> Result goes to RHS    !
        //  diffusion implicit in y and explicit in x      !
        // ================================================!
        // Calculate RHS for Thomas batch algorithm
        prepareRHSjDiff<<<fullDomainGridSize, fullDomainBlockSize>>>(Nx, Ny, phi_x, diff_coeff, tau, d_Vm, d_RHS, d_partRHS);
        CUDA_CALL(cudaDeviceSynchronize());

        // Solve the linear systems
        parallelThomasVertical<<<gridSize_x, THOMAS_KERNEL_BLOCK_SIZE>>>(Nx, Ny, d_RHS, d_la_y, d_lb_y, d_lc_y);
        // parallelThomasVertical<<<gridSize_x, THOMAS_KERNEL_BLOCK_SIZE>>>(Nx, Ny, d_RHS, d_la_y, d_c_prime_y, d_denominator_y);
        CUDA_CALL(cudaDeviceSynchronize());

        // ================================================!
        //  Calculate Vm at n+1 -> Result goes to Vm       !
        //  diffusion implicit in x and explicit in y      !
        // ================================================!
        // Calculate RHS for Thomas batch algorithm
        prepareRHSiDiff<<<fullDomainGridSize, fullDomainBlockSize>>>(Nx, Ny, phi_y, diff_coeff, tau, d_RHS, d_Vm, d_partRHS);
        CUDA_CALL(cudaDeviceSynchronize());

        // Solve the linear systems
        parallelThomasHorizontal<<<gridSize_y, THOMAS_KERNEL_BLOCK_SIZE>>>(Ny, Nx, d_Vm, d_la_x, d_lb_x, d_lc_x);
        // parallelThomasHorizontal<<<gridSize_y, THOMAS_KERNEL_BLOCK_SIZE>>>(Ny, Nx, d_RHS, d_la_x, d_c_prime_x, d_denominator_x);
        CUDA_CALL(cudaDeviceSynchronize());

#endif // SSIADI || THETASSIADI

#ifdef OSADI

        // ================================================!
        //  Calculate Vm at n+1/2 -> Result goes to RHS    !
        // ================================================!
        // Calculate RHS for Thomas batch algorithm
        prepareRHS<<<fullDomainGridSize, fullDomainBlockSize>>>(Nx, Ny, d_Vm, d_partRHS);
        CUDA_CALL(cudaDeviceSynchronize());

        // Solve the linear systems
        parallelThomasVertical<<<gridSize_x, THOMAS_KERNEL_BLOCK_SIZE>>>(Nx, Ny, d_Vm, d_la_y, d_lb_y, d_lc_y);
        CUDA_CALL(cudaDeviceSynchronize());

        // ================================================!
        //  Calculate Vm at n+1 -> Result goes to Vm       !
        // ================================================!
        // Calculate RHS for Thomas batch algorithm
        prepareRHS<<<fullDomainGridSize, fullDomainBlockSize>>>(Nx, Ny, d_Vm, d_partRHS);
        CUDA_CALL(cudaDeviceSynchronize());

        // Solve the linear systems
        parallelThomasHorizontal<<<gridSize_y, THOMAS_KERNEL_BLOCK_SIZE>>>(Ny, Nx, d_Vm, d_la_x, d_lb_x, d_lc_x);
        CUDA_CALL(cudaDeviceSynchronize());

#endif // OSADI

#ifdef FE

        CUDA_CALL(cudaMemcpy(d_Vm, d_partRHS, Nx * Ny * sizeof(real), cudaMemcpyDeviceToDevice));

#endif // FE

#ifdef SAVE_FRAMES

        startSaveFramesTime = omp_get_wtime();

        // If save frames is true and time step is multiple of frame save rate
        if (timeStepCounter % frameSaveRate == 0)
        {
            // Copy memory of d_Vm from device to host Vm
            CUDA_CALL(cudaMemcpy(Vm, d_Vm, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));

            // Save frame
            fprintf(fpFrames, "%lf\n", actualTime);
            saveFrame(fpFrames, Vm, Nx, Ny);
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
            real point_potential = 0.0;
            if (!aux_stim_velocity_flag)
            {
                // Copy memory of d_Vm[begin_point_index] from device to host
                CUDA_CALL(cudaMemcpy(&point_potential, &d_Vm[begin_point_index], sizeof(real), cudaMemcpyDeviceToHost));
#ifdef MV
                point_potential = rescaleVm(point_potential);
#endif // MV
                if (point_potential > 10.0)
                {
                    begin_point_time = actualTime;
                    aux_stim_velocity_flag = true;
                }
            }
            else
            {
                // Copy memory of d_Vm[end_point_index] from device to host
                CUDA_CALL(cudaMemcpy(&point_potential, &d_Vm[end_point_index], sizeof(real), cudaMemcpyDeviceToHost));
#ifdef MV
                point_potential = rescaleVm(point_potential);
#endif // MV
                if (point_potential > 10.0)
                {
                    end_point_time = actualTime;
                    stim_velocity = (end_point - begin_point) / (end_point_time - begin_point_time); // cm/ms
                    stim_velocity = stim_velocity * 10.0;                                            // m/s
                    stim_velocity_measured = true;
                    INFOMSG("Stim velocity (measured from %.2f to %.2f cm) is %.5g m/s\n", begin_point, end_point, stim_velocity);
                }
            }
        }

        finishMeasureVelocityTime = omp_get_wtime();
        elapsedMeasureVelocityTime += finishMeasureVelocityTime - startMeasureVelocityTime;

#endif // MEASURE_VELOCITY

        // Update time step counter
        timeStepCounter++;
    }

#endif // MONODOMAIN

    finishExecutionTime = omp_get_wtime();
    elapsedExecutionTime += finishExecutionTime - startExecutionTime;

    // Copy memory of d_Vm from device to host Vm
    CUDA_CALL(cudaMemcpy(Vm, d_Vm, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));

#ifdef SAVE_FRAMES

    fprintf(fpFrames, "%lf\n", actualTime);
    saveFrame(fpFrames, Vm, Nx, Ny);
    SUCCESSMSG("Frame at time %.2f ms saved to %s\n", actualTime, framesPath);
    fclose(fpFrames);

#endif // SAVE_FRAMES

    printf("Simulation done!\n");
    printf("\n");

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
    fprintf(fpInfos, "BLOCK SIZE = %d x %d threads (%d threads per block)\n", fullDomainBlockSize.x, fullDomainBlockSize.y, fullDomainBlockSize.x * fullDomainBlockSize.y);
    fprintf(fpInfos, "GRID SIZE = %d x %d blocks (total %d blocks, total %d threads)\n", fullDomainGridSize.x, fullDomainGridSize.y, fullDomainGridSize.x * fullDomainGridSize.y, fullDomainGridSize.x * fullDomainGridSize.y * fullDomainBlockSize.x * fullDomainBlockSize.y);

#if defined(SSIADI) || defined(THETASSIADI) || defined(OSADI)

    fprintf(fpInfos, "\n");
    fprintf(fpInfos, "FOR THOMAS KERNEL:\n");
    fprintf(fpInfos, "GRID SIZE FOR X = %d blocks (%d threads per block, total %d threads)\n", gridSize_x, THOMAS_KERNEL_BLOCK_SIZE, gridSize_x * THOMAS_KERNEL_BLOCK_SIZE);
    fprintf(fpInfos, "GRID SIZE FOR Y = %d blocks (%d threads per block, total %d threads)\n", gridSize_y, THOMAS_KERNEL_BLOCK_SIZE, gridSize_y * THOMAS_KERNEL_BLOCK_SIZE);

#endif // SSIADI || THETASSIADI || OSADI

#ifdef MEASURE_VELOCITY

    fprintf(fpInfos, "\n");
    fprintf(fpInfos, "STIMULUS VELOCITY = %.5g m/s\n", stim_velocity);

#endif // MEASURE_VELOCITY

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

    saveFrame(fpLast, Vm, Nx, Ny);

    SUCCESSMSG("Last frame saved to %s\n", lastFrameFilePath);
    fclose(fpLast);

#endif // SAVE_LAST_FRAME

#ifdef SAVE_LAST_STATE

#ifdef AFHN

    char lastFrameFilePathVm[MAX_STRING_SIZE], lastFrameFilePathW[MAX_STRING_SIZE];
    snprintf(lastFrameFilePathVm, MAX_STRING_SIZE * sizeof(char), "%s/lastframeVm.txt", pathToSaveData);
    snprintf(lastFrameFilePathW, MAX_STRING_SIZE * sizeof(char), "%s/lastframeW.txt", pathToSaveData);

    CUDA_CALL(cudaMemcpy(W, d_W, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));

    FILE *fpLastVm = fopen(lastFrameFilePathVm, "w");
    FILE *fpLastW = fopen(lastFrameFilePathW, "w");

    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            int index = i * Nx + j;
            fprintf(fpLastVm, "%e ", Vm[index]);
            fprintf(fpLastW, "%e ", W[index]);
        }
        fprintf(fpLastVm, "\n");
        fprintf(fpLastW, "\n");
    }

    SUCCESSMSG("Last Vm frame saved to %s\n", lastFrameFilePathVm);
    SUCCESSMSG("Last W frame saved to %s\n", lastFrameFilePathW);
    fclose(fpLastVm);
    fclose(fpLastW);

#endif // AFHN

#ifdef TT2

    char lastFrameFilePathVm[MAX_STRING_SIZE], lastFrameFilePathX_r1[MAX_STRING_SIZE], lastFrameFilePathX_r2[MAX_STRING_SIZE], lastFrameFilePathX_s[MAX_STRING_SIZE],
        lastFrameFilePathm[MAX_STRING_SIZE], lastFrameFilePathh[MAX_STRING_SIZE], lastFrameFilePathj[MAX_STRING_SIZE], lastFrameFilePathd[MAX_STRING_SIZE], lastFrameFilePathf[MAX_STRING_SIZE],
        lastFrameFilePathf2[MAX_STRING_SIZE], lastFrameFilePathfCaSS[MAX_STRING_SIZE], lastFrameFilePaths[MAX_STRING_SIZE], lastFrameFilePathr[MAX_STRING_SIZE],
        lastFrameFilePathR_prime[MAX_STRING_SIZE], lastFrameFilePathCa_i[MAX_STRING_SIZE], lastFrameFilePathCa_SR[MAX_STRING_SIZE], lastFrameFilePathCa_SS[MAX_STRING_SIZE],
        lastFrameFilePathNa_i[MAX_STRING_SIZE], lastFrameFilePathK_i[MAX_STRING_SIZE];
    snprintf(lastFrameFilePathVm, MAX_STRING_SIZE * sizeof(char), "%s/lastframeVm.txt", pathToSaveData);
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

    FILE *fpLastVm = fopen(lastFrameFilePathVm, "w");
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
            fprintf(fpLastVm, "%e ", Vm[index]);
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
        fprintf(fpLastVm, "\n");
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

    SUCCESSMSG("Last Vm frame saved to %s\n", lastFrameFilePathVm);
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

    fclose(fpLastVm);
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

#ifdef MV

    char lastFrameFilePathVm[MAX_STRING_SIZE], lastFrameFilePathv[MAX_STRING_SIZE], lastFrameFilePathw[MAX_STRING_SIZE], lastFrameFilePaths[MAX_STRING_SIZE];
    snprintf(lastFrameFilePathVm, MAX_STRING_SIZE * sizeof(char), "%s/lastframeVm.txt", pathToSaveData);
    snprintf(lastFrameFilePathv, MAX_STRING_SIZE * sizeof(char), "%s/lastframev.txt", pathToSaveData);
    snprintf(lastFrameFilePathw, MAX_STRING_SIZE * sizeof(char), "%s/lastframew.txt", pathToSaveData);
    snprintf(lastFrameFilePaths, MAX_STRING_SIZE * sizeof(char), "%s/lastframes.txt", pathToSaveData);

    CUDA_CALL(cudaMemcpy(v, d_v, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(w, d_w, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(s, d_s, Nx * Ny * sizeof(real), cudaMemcpyDeviceToHost));

    FILE *fpLastVm = fopen(lastFrameFilePathVm, "w");
    FILE *fpLastv = fopen(lastFrameFilePathv, "w");
    FILE *fpLastw = fopen(lastFrameFilePathw, "w");
    FILE *fpLasts = fopen(lastFrameFilePaths, "w");

    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            int index = i * Nx + j;
            fprintf(fpLastVm, "%e ", Vm[index]);
            fprintf(fpLastv, "%e ", v[index]);
            fprintf(fpLastw, "%e ", w[index]);
            fprintf(fpLasts, "%e ", s[index]);
        }
        fprintf(fpLastVm, "\n");
        fprintf(fpLastv, "\n");
        fprintf(fpLastw, "\n");
        fprintf(fpLasts, "\n");
    }

    SUCCESSMSG("Last Vm frame saved to %s\n", lastFrameFilePathVm);
    SUCCESSMSG("Last v frame saved to %s\n", lastFrameFilePathv);
    SUCCESSMSG("Last w frame saved to %s\n", lastFrameFilePathw);
    SUCCESSMSG("Last s frame saved to %s\n", lastFrameFilePaths);

#endif // MV

#endif // SAVE_LAST_STATE

    // Free memory
    free(time);
    free(pathToSaveData);

    CUDA_CALL(cudaFree(d_partRHS));

#if defined(SSIADI) || defined(THETASSIADI) || defined(OSADI)

    CUDA_CALL(cudaFree(d_c_prime_x));
    CUDA_CALL(cudaFree(d_denominator_x));
    CUDA_CALL(cudaFree(d_c_prime_y));
    CUDA_CALL(cudaFree(d_denominator_y));

    CUDA_CALL(cudaFree(d_la_x));
    CUDA_CALL(cudaFree(d_lb_x));
    CUDA_CALL(cudaFree(d_lc_x));
    CUDA_CALL(cudaFree(d_la_y));
    CUDA_CALL(cudaFree(d_lb_y));
    CUDA_CALL(cudaFree(d_lc_y));

#ifndef OSADI

    CUDA_CALL(cudaFree(d_RHS));

#endif // OSADI

#endif // SSIADI || THETASSIADI || OSADI

#ifdef MONODOMAIN

    free(stimuli);
    CUDA_CALL(cudaFree(d_stimuli));

#ifdef AFHN

    free(Vm);
    free(W);
    CUDA_CALL(cudaFree(d_Vm));
    CUDA_CALL(cudaFree(d_W));

#endif // AFHN

#ifdef TT2

    free(Vm);
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
    CUDA_CALL(cudaFree(d_Vm));
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

#ifdef MV

    free(Vm);
    free(v);
    free(w);
    free(s);
    CUDA_CALL(cudaFree(d_Vm));
    CUDA_CALL(cudaFree(d_v));
    CUDA_CALL(cudaFree(d_w));
    CUDA_CALL(cudaFree(d_s));

#endif // MV
#endif // MONODOMAIN

    // Reset device
    CUDA_CALL(cudaDeviceReset());

    return;
}

#endif // GPU_METHODS_H