#ifndef CPU_METHODS_H
#define CPU_METHODS_H

#include "auxfuncs.h"

void runSimulation(char *method, real delta_t, real delta_x, real theta)
{
    // Number of steps
    int N = round(L / delta_x) + 1;     // Spatial steps (square tissue)
    int M = round(totalTime / delta_t); // Number of time steps

    // Allocate and populate time array
    real *time = (real *)malloc(M * sizeof(real));
    initializeTimeArray(time, M, delta_t);
#ifndef CABLEEQ
    // Allocate 2D arrays for variables
    real **V, **Vtilde, **RHS, **partRHS, **exact;
    V = (real **)malloc(N * sizeof(real *));
    Vtilde = (real **)malloc(N * sizeof(real *));
    RHS = (real **)malloc(N * sizeof(real *));
    partRHS = (real **)malloc(N * sizeof(real *));
    exact = (real **)malloc(N * sizeof(real *));
#else // if def CABLEEQ
    // Allocate 1D arrays for variables
    real *V, *Vtilde, *RHS, *partRHS, *exact, *AP;
    V = (real *)malloc(N * sizeof(real));
    Vtilde = (real *)malloc(N * sizeof(real));
    RHS = (real *)malloc(N * sizeof(real));
    partRHS = (real *)malloc(N * sizeof(real));
    exact = (real *)malloc(N * sizeof(real));
    AP = (real *)malloc(M * sizeof(real));
#endif // not CABLEEQ
    // Aux variables for the Linear System resolution (Thomas)
    real *c_prime = (real *)malloc(N * sizeof(real)); // aux for Thomas
    real *d_prime = (real *)malloc(N * sizeof(real)); // aux for Thomas
    real *LS_b = (real *)malloc(N * sizeof(real));
    real *result = (real *)malloc(N * sizeof(real));
#ifdef MONODOMAIN
#ifdef AFHN
    real **W = (real **)malloc(N * sizeof(real *));
#endif // AFHN
#ifdef TT2
    real **X_r1, **X_r2, **X_s, **m, **h, **_j, **d, **f, **f2, **fCaSS, **s, **r, **Ca_i, **Ca_SR, **Ca_SS, **R_prime, **Na_i, **K_i;
    X_r1 = (real **)malloc(N * sizeof(real *));
    X_r2 = (real **)malloc(N * sizeof(real *));
    X_s = (real **)malloc(N * sizeof(real *));
    m = (real **)malloc(N * sizeof(real *));
    h = (real **)malloc(N * sizeof(real *));
    _j = (real **)malloc(N * sizeof(real *));
    d = (real **)malloc(N * sizeof(real *));
    f = (real **)malloc(N * sizeof(real *));
    f2 = (real **)malloc(N * sizeof(real *));
    fCaSS = (real **)malloc(N * sizeof(real *));
    s = (real **)malloc(N * sizeof(real *));
    r = (real **)malloc(N * sizeof(real *));
    Ca_i = (real **)malloc(N * sizeof(real *));
    Ca_SR = (real **)malloc(N * sizeof(real *));
    Ca_SS = (real **)malloc(N * sizeof(real *));
    R_prime = (real **)malloc(N * sizeof(real *));
    Na_i = (real **)malloc(N * sizeof(real *));
    K_i = (real **)malloc(N * sizeof(real *));
#endif // TT2
#endif // MONODOMAIN
#ifdef CABLEEQ
#ifdef AFHN
    real *W = (real *)malloc(N * sizeof(real));
#endif // AFHN
#ifdef TT2
    real *X_r1, *X_r2, *X_s, *m, *h, *_j, *d, *f, *f2, *fCaSS, *s, *r, *Ca_i, *Ca_SR, *Ca_SS, *R_prime, *Na_i, *K_i;
    X_r1 = (real *)malloc(N * sizeof(real));
    X_r2 = (real *)malloc(N * sizeof(real));
    X_s = (real *)malloc(N * sizeof(real));
    m = (real *)malloc(N * sizeof(real));
    h = (real *)malloc(N * sizeof(real));
    _j = (real *)malloc(N * sizeof(real));
    d = (real *)malloc(N * sizeof(real));
    f = (real *)malloc(N * sizeof(real));
    f2 = (real *)malloc(N * sizeof(real));
    fCaSS = (real *)malloc(N * sizeof(real));
    s = (real *)malloc(N * sizeof(real));
    r = (real *)malloc(N * sizeof(real));
    Ca_i = (real *)malloc(N * sizeof(real));
    Ca_SR = (real *)malloc(N * sizeof(real));
    Ca_SS = (real *)malloc(N * sizeof(real));
    R_prime = (real *)malloc(N * sizeof(real));
    Na_i = (real *)malloc(N * sizeof(real));
    K_i = (real *)malloc(N * sizeof(real));
#endif // TT2
#endif // CABLEEQ
#ifndef CABLEEQ
    for (int i = 0; i < N; i++)
    {
        V[i] = (real *)malloc(N * sizeof(real));
        Vtilde[i] = (real *)malloc(N * sizeof(real));
        RHS[i] = (real *)malloc(N * sizeof(real));
        partRHS[i] = (real *)malloc(N * sizeof(real));
        exact[i] = (real *)malloc(N * sizeof(real));
#ifdef MONODOMAIN
#ifdef AFHN
        W[i] = (real *)malloc(N * sizeof(real));
#endif // AFHN
#ifdef TT2
        X_r1[i] = (real *)malloc(N * sizeof(real));
        X_r2[i] = (real *)malloc(N * sizeof(real));
        X_s[i] = (real *)malloc(N * sizeof(real));
        m[i] = (real *)malloc(N * sizeof(real));
        h[i] = (real *)malloc(N * sizeof(real));
        _j[i] = (real *)malloc(N * sizeof(real));
        d[i] = (real *)malloc(N * sizeof(real));
        f[i] = (real *)malloc(N * sizeof(real));
        f2[i] = (real *)malloc(N * sizeof(real));
        fCaSS[i] = (real *)malloc(N * sizeof(real));
        s[i] = (real *)malloc(N * sizeof(real));
        r[i] = (real *)malloc(N * sizeof(real));
        Ca_i[i] = (real *)malloc(N * sizeof(real));
        Ca_SR[i] = (real *)malloc(N * sizeof(real));
        Ca_SS[i] = (real *)malloc(N * sizeof(real));
        R_prime[i] = (real *)malloc(N * sizeof(real));
        Na_i[i] = (real *)malloc(N * sizeof(real));
        K_i[i] = (real *)malloc(N * sizeof(real));
#endif // TT2
#endif // MONODOMAIN
    }
#endif // not CABLEEQ

#ifdef CONVERGENCE_ANALYSIS
    initialize2DVariableWithExactSolution(V, N, delta_x);
#else // if not def CONVERGENCE_ANALYSIS
#ifndef CABLEEQ
    initialize2DVariableWithValue(V, N, V_init);
#else // if def CABLEEQ
    initialize1DVariableWithValue(V, N, V_init);
#endif // not CABLEEQ
#endif // CONVERGENCE_ANALYSIS

#ifdef MONODOMAIN
#ifdef AFHN
    initialize2DVariableWithValue(W, N, W_init);
#endif // AFHN
#ifdef TT2
    initialize2DVariableWithValue(X_r1, N, X_r1_init);
    initialize2DVariableWithValue(X_r2, N, X_r2_init);
    initialize2DVariableWithValue(X_s, N, X_s_init);
    initialize2DVariableWithValue(m, N, m_init);
    initialize2DVariableWithValue(h, N, h_init);
    initialize2DVariableWithValue(_j, N, j_init);
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
#ifdef CABLEEQ
#ifdef AFHN
    initialize1DVariableWithValue(W, N, W_init);
#endif // AFHN
#ifdef TT2
    initialize1DVariableWithValue(X_r1, N, X_r1_init);
    initialize1DVariableWithValue(X_r2, N, X_r2_init);
    initialize1DVariableWithValue(X_s, N, X_s_init);
    initialize1DVariableWithValue(m, N, m_init);
    initialize1DVariableWithValue(h, N, h_init);
    initialize1DVariableWithValue(_j, N, j_init);
    initialize1DVariableWithValue(d, N, d_init);
    initialize1DVariableWithValue(f, N, f_init);
    initialize1DVariableWithValue(f2, N, f2_init);
    initialize1DVariableWithValue(fCaSS, N, fCaSS_init);
    initialize1DVariableWithValue(s, N, s_init);
    initialize1DVariableWithValue(r, N, r_init);
    initialize1DVariableWithValue(Ca_i, N, Ca_i_init);
    initialize1DVariableWithValue(Ca_SR, N, Ca_SR_init);
    initialize1DVariableWithValue(Ca_SS, N, Ca_SS_init);
    initialize1DVariableWithValue(R_prime, N, R_prime_init);
    initialize1DVariableWithValue(Na_i, N, Na_i_init);
    initialize1DVariableWithValue(K_i, N, K_i_init);
#endif // TT2
#endif // CABLEEQ

#ifdef INIT_WITH_SPIRAL
    char* reference_dt = "0.00010";
    char* reference_dx = "0.00050";
    real real_ref_dx = 0.0005f;
    char *pathToSpiralFiles = (char *)malloc(MAX_STRING_SIZE * sizeof(char));
    snprintf(pathToSpiralFiles, MAX_STRING_SIZE * sizeof(char), "./spiral_files/%s/%s/%s/lastV_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize2DVariableFromFile(V, N, pathToSpiralFiles, delta_x, "V", real_ref_dx);
    snprintf(pathToSpiralFiles, MAX_STRING_SIZE * sizeof(char), "./spiral_files/%s/%s/%s/lastW_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize2DVariableFromFile(W, N, pathToSpiralFiles, delta_x, "W", real_ref_dx);
    free(pathToSpiralFiles);
#endif // INIT_WITH_SPIRAL

#ifdef RESTORE_STATE_AND_SHIFT
    // Initialize variables with a solution
    char* reference_dt = "0.00010";
    char* reference_dx = "0.00050";
    real real_ref_dx = 0.0005f;
    char *pathToRestoreStateFiles = (char *)malloc(MAX_STRING_SIZE * sizeof(char));

    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastV_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    #ifdef CABLEEQ
    initialize1DVariableFromFile(V, N, pathToRestoreStateFiles, delta_x, "V", real_ref_dx);
    #ifdef AFHN
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastW_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize1DVariableFromFile(W, N, pathToRestoreStateFiles, delta_x, "W", real_ref_dx);
    #endif // AFHN
    #ifdef TT2
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastX_r1_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize1DVariableFromFile(X_r1, N, pathToRestoreStateFiles, delta_x, "X_r1", real_ref_dx);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastX_r2_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize1DVariableFromFile(X_r2, N, pathToRestoreStateFiles, delta_x, "X_r2", real_ref_dx);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastX_s_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize1DVariableFromFile(X_s, N, pathToRestoreStateFiles, delta_x, "X_s", real_ref_dx);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastm_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize1DVariableFromFile(m, N, pathToRestoreStateFiles, delta_x, "m", real_ref_dx);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lasth_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize1DVariableFromFile(h, N, pathToRestoreStateFiles, delta_x, "h", real_ref_dx);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastj_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize1DVariableFromFile(_j, N, pathToRestoreStateFiles, delta_x, "j", real_ref_dx);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastd_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize1DVariableFromFile(d, N, pathToRestoreStateFiles, delta_x, "d", real_ref_dx);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastf_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize1DVariableFromFile(f, N, pathToRestoreStateFiles, delta_x, "f", real_ref_dx);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastf2_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize1DVariableFromFile(f2, N, pathToRestoreStateFiles, delta_x, "f2", real_ref_dx);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastfCaSS_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize1DVariableFromFile(fCaSS, N, pathToRestoreStateFiles, delta_x, "fCaSS", real_ref_dx);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lasts_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize1DVariableFromFile(s, N, pathToRestoreStateFiles, delta_x, "s", real_ref_dx);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastr_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize1DVariableFromFile(r, N, pathToRestoreStateFiles, delta_x, "r", real_ref_dx);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastCa_i_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize1DVariableFromFile(Ca_i, N, pathToRestoreStateFiles, delta_x, "Ca_i", real_ref_dx);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastCa_SR_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize1DVariableFromFile(Ca_SR, N, pathToRestoreStateFiles, delta_x, "Ca_SR", real_ref_dx);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastCa_SS_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize1DVariableFromFile(Ca_SS, N, pathToRestoreStateFiles, delta_x, "Ca_SS", real_ref_dx);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastR_prime_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize1DVariableFromFile(R_prime, N, pathToRestoreStateFiles, delta_x, "R_prime", real_ref_dx);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastNa_i_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize1DVariableFromFile(Na_i, N, pathToRestoreStateFiles, delta_x, "Na_i", real_ref_dx);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastK_i_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize1DVariableFromFile(K_i, N, pathToRestoreStateFiles, delta_x, "K_i", real_ref_dx);
    #endif // TT2
    #endif // CABLEEQ
    #ifdef MONODOMAIN
    initialize2DVariableFromFile(V, N, pathToRestoreStateFiles, delta_x, "V", real_ref_dx);
    #ifdef AFHN
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE * sizeof(char), "./restore_state/%s/%s/%s/lastW_%s_%s.txt", REAL_TYPE, PROBLEM, CELL_MODEL, reference_dt, reference_dx);
    initialize2DVariableFromFile(W, N, pathToRestoreStateFiles, delta_x, "W", real_ref_dx);
    #endif // AFHN
    #endif MONODOMAIN
    free(pathToRestoreStateFiles);

    // Shift variables
    real lengthToShift = 0.1f;

    #ifdef CABLEEQ
    shift1DVariableToLeft(V, N, lengthToShift, delta_x, V_init, "V");
    #ifdef AFHN
    shift1DVariableToLeft(W, N, lengthToShift, delta_x, W_init, "W");
    #endif // AFHN
    #ifdef TT2
    shift1DVariableToLeft(X_r1, N, lengthToShift, delta_x, X_r1_init, "X_r1");
    shift1DVariableToLeft(X_r2, N, lengthToShift, delta_x, X_r2_init, "X_r2");
    shift1DVariableToLeft(X_s, N, lengthToShift, delta_x, X_s_init, "X_s");
    shift1DVariableToLeft(m, N, lengthToShift, delta_x, m_init, "m");
    shift1DVariableToLeft(h, N, lengthToShift, delta_x, h_init, "h");
    shift1DVariableToLeft(_j, N, lengthToShift, delta_x, j_init, "j");
    shift1DVariableToLeft(d, N, lengthToShift, delta_x, d_init, "d");
    shift1DVariableToLeft(f, N, lengthToShift, delta_x, f_init, "f");
    shift1DVariableToLeft(f2, N, lengthToShift, delta_x, f2_init, "f2");
    shift1DVariableToLeft(fCaSS, N, lengthToShift, delta_x, fCaSS_init, "fCaSS");
    shift1DVariableToLeft(s, N, lengthToShift, delta_x, s_init, "s");
    shift1DVariableToLeft(r, N, lengthToShift, delta_x, r_init, "r");
    shift1DVariableToLeft(Ca_i, N, lengthToShift, delta_x, Ca_i_init, "Ca_i");
    shift1DVariableToLeft(Ca_SR, N, lengthToShift, delta_x, Ca_SR_init, "Ca_SR");
    shift1DVariableToLeft(Ca_SS, N, lengthToShift, delta_x, Ca_SS_init, "Ca_SS");
    shift1DVariableToLeft(R_prime, N, lengthToShift, delta_x, R_prime_init, "R_prime");
    shift1DVariableToLeft(Na_i, N, lengthToShift, delta_x, Na_i_init, "Na_i");
    shift1DVariableToLeft(K_i, N, lengthToShift, delta_x, K_i_init, "K_i");
    #endif // TT2
    #endif // CABLEEQ
    #ifdef MONODOMAIN
    shift2DVariableToLeft(V, N, lengthToShift, delta_x, V_init, "V");
    #ifdef AFHN
    shift2DVariableToLeft(W, N, lengthToShift, delta_x, W_init, "W");
    #endif // AFHN
    #endif // MONODOMAIN
#endif // RESTORE_STATE_AND_SHIFT

    // Auxiliary arrays for Thomas algorithm
    real *la = (real *)malloc(N * sizeof(real)); // subdiagonal
    real *lb = (real *)malloc(N * sizeof(real)); // diagonal
    real *lc = (real *)malloc(N * sizeof(real)); // superdiagonal

    // Populate auxiliary arrays for Thomas algorithm
    real phi = delta_t / (delta_x * delta_x);
    #ifndef TT2
    real diff_coeff = sigma / (Cm * chi);
    #else
    real diff_coeff = sigma / chi;
    #endif // not TT2
    if (strcmp(method, "ADI") == 0 || strcmp(method, "SSI-ADI") == 0)
    {
        populateDiagonalThomasAlgorithm(la, lb, lc, N, 0.5f*phi*diff_coeff);
    }
    else if (strstr(method, "theta") != NULL)
    {
        populateDiagonalThomasAlgorithm(la, lb, lc, N, theta*phi*diff_coeff);
    }

#ifndef CONVERGENCE_ANALYSIS
#if defined(MONODOMAIN) || defined(CABLEEQ)
    // Allocate array for the stimuli
    Stimulus *stimuli = (Stimulus *)malloc(numberOfStimuli * sizeof(Stimulus));
    populateStimuli(stimuli, delta_x);
#endif // MONODOMAIN || CABLEEQ
#endif // not CONVERGENCE_ANALYSIS

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

#ifdef CABLEEQ
    // Choose cell at 0.5 cm to measure Action Potential
    int APCellIndex = round(0.5f/delta_x);
#endif // CABLEEQ

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
#ifndef CABLEEQ
    if (strcmp(method, "ADI") == 0)
    {
        while (timeStepCounter < M)
        {
            // Get time step
            actualTime = time[timeStepCounter];

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
#if defined(LINMONO) || defined(DIFF)
                    LS_b[i] = 0.5f * phi * V[i][lim(j - 1, N)] + (1.0f - 2.0f * 0.5f * phi) * V[i][j] + 0.5f * phi * V[i][lim(j + 1, N)] + 0.5f * delta_t * forcingTerm(x, y, actualTime);
#endif // LINMONO || DIFF
                }

                tridiag(la, lb, lc, c_prime, d_prime, N, LS_b, result);
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
#if defined(LINMONO) || defined(DIFF)
                    LS_b[j] = 0.5f * phi * RHS[lim(i - 1, N)][j] + (1.0f - 2.0f * 0.5f * phi) * RHS[i][j] + 0.5f * phi * RHS[lim(i + 1, N)][j] + 0.5f * delta_t * forcingTerm(x, y, actualTime + delta_t);
#endif // LINMONO || DIFF
                }

                tridiag(la, lb, lc, c_prime, d_prime, N, LS_b, result);
                for (int j = 0; j < N; j++)
                {
                    V[i][j] = result[j];
                }
            }

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
            //  Calculate Approxs. and Update ODEs             !
            // ================================================!
            real x, y;
            real diff_term = 0.0f;
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    x = j * delta_x;
                    y = i * delta_x;

#ifdef LINMONO
                    diff_term = diff_coeff * 0.5f * phi * (V[i][lim(j - 1, N)] + V[lim(i - 1, N)][j] - 4.0f * V[i][j] + V[i][lim(j + 1, N)] + V[lim(i + 1, N)][j]);
                    real for_term = forcingTerm(x, y, actualTime + (0.5f * delta_t)) / (chi * Cm);
                    real reac_term = G * V[i][j] / Cm;
                    real actualVtilde = V[i][j] + diff_term + (0.5f * delta_t * (for_term - reac_term));

                    // Preparing part of the RHS of the following linear systems
                    real reac_tilde_term = G * actualVtilde / Cm;
                    partRHS[i][j] = delta_t * (for_term - reac_tilde_term);
#endif // LINMONO

#ifdef MONODOMAIN
#ifdef AFHN
                    real actualV = V[i][j];
                    real actualW = W[i][j];
                    diff_term = diff_coeff * phi * (V[i][lim(j - 1, N)] + V[lim(i - 1, N)][j] - 4.0f * actualV + V[i][lim(j + 1, N)] + V[lim(i + 1, N)][j]);
                    real RHS_V_term = RHS_V(actualV, actualW) / (Cm * chi);
#ifdef CONVERGENCE_ANALYSIS
                    real for_term = forcingTerm(x, y, actualTime + (0.5f * delta_t), actualW) / (chi * Cm);
#endif // CONVERGENCE_ANALYSIS

#ifndef CONVERGENCE_ANALYSIS
                    // Stimulation
                    real stim = 0.0f;
                    for (int si = 0; si < numberOfStimuli; si++)
                    {
                        if (actualTime >= stimuli[si].begin && actualTime <= stimuli[si].begin + stimuli[si].duration && j >= stimuli[si].xMinDisc && j <= stimuli[si].xMaxDisc && i >= stimuli[si].yMinDisc && i <= stimuli[si].yMaxDisc)
                        {
                            stim = stimuli[si].strength;
                            break;
                        }
                    }

                    real actualVtilde = actualV + 0.5f * diff_term + (0.5f * delta_t * (stim - RHS_V_term));
#else
                    real actualVtilde = actualV + 0.5f * diff_term + (0.5f * delta_t * (for_term - RHS_V_term));
#endif // not CONVERGENCE_ANALYSIS
                    // Calculate approximation for state variables
                    real Wtilde = actualW + (0.5f * delta_t * RHS_W(actualV, actualW));

                    // Preparing part of the RHS of the following linear systems
                    real RHS_Vtilde_term = RHS_V(actualVtilde, Wtilde) / (Cm * chi);

                    #ifndef CONVERGENCE_ANALYSIS
                    partRHS[i][j] = delta_t * (stim - RHS_Vtilde_term);
                    #else
                    partRHS[i][j] = delta_t * (for_term - RHS_Vtilde_term);
                    #endif // not CONVERGENCE_ANALYSIS

                    // Update state variables with RK2 -> Wn+1 = Wn + dt*R(V*, W*)
                    W[i][j] = actualW + delta_t * RHS_W(actualVtilde, Wtilde);
#endif // AFHN
#ifdef TT2
                    real actualV = V[i][j];
                    real actualX_r1 = X_r1[i][j];
                    real actualX_r2 = X_r2[i][j];
                    real actualX_s = X_s[i][j];
                    real actualm = m[i][j];
                    real actualh = h[i][j];
                    real actualj = _j[i][j];
                    real actuald = d[i][j];
                    real actualf = f[i][j];
                    real actualf2 = f2[i][j];
                    real actualfCaSS = fCaSS[i][j];
                    real actuals = s[i][j];
                    real actualr = r[i][j];
                    real actualCa_i = Ca_i[i][j];
                    real actualCa_SR = Ca_SR[i][j];
                    real actualCa_SS = Ca_SS[i][j];
                    real actualR_prime = R_prime[i][j];
                    real actualNa_i = Na_i[i][j];
                    real actualK_i = K_i[i][j];

                    diff_term = diff_coeff * phi * (V[i][lim(j - 1, N)] + V[lim(i - 1, N)][j] - 4.0f * actualV + V[i][lim(j + 1, N)] + V[lim(i + 1, N)][j]);

                    // Auxiliary variables for the calculation of the RHS of the ODEs
                    real VmENa = actualV - (RTONF * log(Na_o / actualNa_i));
                    real E_K = (RTONF * log(K_o / actualK_i));
                    real VmEK = actualV - E_K;
                    real alpha_K1 = 0.1f / (1.0f + exp(0.06f * (actualV - E_K - 200.0f)));
                    real beta_K1 = (3.0f * exp(0.0002f * (actualV - E_K + 100.0f)) + exp(0.1f * (actualV - E_K - 10.0f))) / (1.0f + exp(-0.5f * (actualV - E_K)));
                    real E_Ks = (RTONF * log((K_o + p_KNa * Na_o) / (actualK_i + p_KNa * actualNa_i)));
                    real E_Ca = 0.5f * RTONF * log(Ca_o / actualCa_i);

                    // Currents
                    real INa = G_Na * (actualm * actualm * actualm) * actualh * actualj * VmENa;
                    real IbNa = G_bNa * VmENa;
                    real IK1 = G_K1 * (alpha_K1 / (alpha_K1 + beta_K1)) * VmEK;
                    real Ito = G_to * actualr * actuals * VmEK;
                    real IKr = G_Kr * sqrt(K_o / 5.4f) * actualX_r1 * actualX_r2 * VmEK;
                    real IKs = G_Ks * actualX_s * actualX_s * (actualV - E_Ks);
                    real ICaL; // !!!
                    (actualV < 15.0f - 1.0e-5f)
                        ? (ICaL = G_CaL * actuald * actualf * actualf2 * actualfCaSS * 4.0f * (actualV - 15.0f) * (F * F) * (0.25f * actualCa_SS * exp(2.0f * (actualV - 15.0f) * FONRT) - Ca_o) / (R * T * (exp(2.0f * (actualV - 15.0f) * FONRT) - 1.0f)))
                        : (ICaL = G_CaL * actuald * actualf * actualf2 * actualfCaSS * 2.0f * F * (0.25f * actualCa_SS - Ca_o));
                    real INaK = ((((p_KNa * K_o) / (K_o + K_mK)) * actualNa_i) / (actualNa_i + K_mNa)) / (1.0f + (0.1245f * exp(((-0.1f) * actualV * FONRT))) + (0.0353f * exp(((-actualV) * FONRT))));
                    real INaCa; // !!!
                    INaCa = (k_NaCa * ((exp((gamma_I_NaCa * actualV * FONRT)) * (actualNa_i * actualNa_i * actualNa_i) * Ca_o) - (exp(((gamma_I_NaCa - 1.0f) * actualV * FONRT)) * (Na_o * Na_o * Na_o) * actualCa_i * alpha))) / (((K_mNa_i * K_mNa_i * K_mNa_i) + (Na_o * Na_o * Na_o)) * (K_mCa + Ca_o) * (1.0f + (k_sat * exp(((gamma_I_NaCa)*actualV * FONRT)))));
                    real IpCa = (G_pCa * actualCa_i) / (K_pCa + actualCa_i);
                    real IpK = (G_pK * VmEK) / (1.0f + exp((25.0f - actualV) / 5.98f));
                    real IbCa = G_bCa * (actualV - E_Ca);

                    // RHS_V at actual time
                    real RHS_V_term = INa + IbNa + IK1 + Ito + IKr + IKs + ICaL + INaK + INaCa + IpCa + IpK + IbCa;

                    // Stimulation
                    real stim = 0.0f;
                    for (int si = 0; si < numberOfStimuli; si++)
                    {
                        if (actualTime >= stimuli[si].begin && actualTime <= stimuli[si].begin + stimuli[si].duration && j >= stimuli[si].xMinDisc && j <= stimuli[si].xMaxDisc && i >= stimuli[si].yMinDisc && i <= stimuli[si].yMaxDisc)
                        {
                            stim = stimuli[si].strength;
                            break;
                        }
                    }

                    // Calculate Vtilde -> utilde = u^n + 0.5 * dt * (A*u^n + R(u^n))
                    real actualVtilde = actualV + 0.5f * diff_term + (0.5f * delta_t * (stim - RHS_V_term));

                    // Preparing part of the RHS of the following linear systems
                    // Calculate approximation for state variables
                    // Rush-Larsen method - auxiliary variables
                    real X_r1_inf = 1.0f / (1.0f + exp((-26.0f - actualV) / 7.0f));
                    real alpha_X_r1 = 450.0f / (1.0f + exp((-45.0f - actualV) / 10.0f));
                    real beta_X_r1 = 6.0f / (1.0f + exp((30.0f + actualV) / 11.5f));

                    real X_r2_inf = 1.0f / (1.0f + exp((actualV + 88.0f) / 24.0f));
                    real alpha_X_r2 = 3.0f / (1.0f + exp((-60.0f - actualV) / 20.0f));
                    real beta_X_r2 = 1.12f / (1.0f + exp((actualV - 60.0f) / 20.0f));

                    real X_s_inf = 1.0f / (1.0f + exp((-5.0f - actualV) / 14.0f));
                    real alpha_X_s = 1400.0f / sqrt(1.0f + exp((5.0f - actualV) / 6.0f));
                    real beta_X_s = 1.0f / (1.0f + exp((-35.0f + actualV) / 15.0f));

                    real m_inf = 1.0f / ((1.0f + exp((-56.86 - actualV) / 9.03f)) * (1.0f + exp((-56.86f - actualV) / 9.03f)));
                    real alpha_m = 1.0f / (1.0f + exp((-60.0f - actualV) / 5.0f));
                    real beta_m = 0.1f / (1.0f + exp((actualV + 35.0f) / 5.0f)) + (0.1f / (1.0f + exp((actualV - 50.0f) / 200.0f)));

                    real h_inf = 1.0f / ((1.0f + exp((actualV + 71.55f) / 7.43f)) * (1.0f + exp((actualV + 71.55f) / 7.43f)));
                    real alpha_h;
                    (actualV < -40.0f)
                        ? (alpha_h = 0.057f * exp(-(actualV + 80.0f) / 6.8f))
                        : (alpha_h = 0.0f);
                    real beta_h;
                    (actualV < -40.0f)
                        ? (beta_h = 2.7f * exp(0.079f * actualV) + 3.1f * 1.0e5f * exp(0.3485f * actualV))
                        : (beta_h = 0.77f / (0.13f * (1.0f + exp((actualV + 10.66f) / -11.1f))));

                    real j_inf = 1.0f / ((1.0f + exp((actualV + 71.55f) / 7.43f)) * (1.0f + exp((actualV + 71.55f) / 7.43f)));
                    real alpha_j;
                    (actualV < -40.0f)
                        ? (alpha_j = ((-25428.0f * exp(0.2444f * actualV) - (6.948e-6f * exp((-0.04391f) * actualV))) * (actualV + 37.78f)) / (1.0f + exp(0.311f * (actualV + 79.23f))))
                        : (alpha_j = 0.0f);
                    real beta_j;
                    (actualV < -40.0f)
                        ? (beta_j = (0.02424f * exp(-0.01052f * actualV)) / (1.0f + exp(-0.1378f * (actualV + 40.14f))))
                        : (beta_j = (0.6f * exp(0.057f * actualV)) / (1.0f + exp(-0.1f * (actualV + 32.0f))));

                    real inf = 1.0f / (1.0f + exp((-8.0f - actualV) / 7.5f));
                    real alpha_d = 1.4f / (1.0f + exp((-35.0f - actualV) / 13.0f)) + 0.25f;
                    real beta_d = 1.4f / (1.0f + exp((actualV + 5.0f) / 5.0f));
                    real gamma_d = 1.0f / (1.0f + exp((50.0f - actualV) / 20.0f));

                    real f_inf = 1.0f / (1.0f + exp((actualV + 20.0f) / 7.0f));
                    real alpha_f = 1102.5f * exp(-(actualV + 27.0f) * (actualV + 27.0f) / 225.0f);
                    real beta_f = 200.0f / (1.0f + exp((13.0f - actualV) / 10.0f));
                    real gamma_f = 180.0f / (1.0f + exp((actualV + 30.0f) / 10.0f)) + 20.0f;

                    real f2_inf = 0.67f / (1.0f + exp((actualV + 35.0f) / 7.0f)) + 0.33f;
                    real alpha_f2; // !!!
                    alpha_f2 = 562.0f * exp(-(actualV + 27.0f) * (actualV + 27.0f) / 240.0f);
                    real beta_f2 = 31.0f / (1.0f + exp((25.0f - actualV) / 10.0f));
                    real gamma_f2; // !!!
                    gamma_f2 = 80.0f / (1.0f + exp((30.0f + actualV) / 10.0f));

                    real fCaSS_inf = 0.6f / (1.0f + (actualCa_SS * actualCa_SS * 400.0f)) + 0.4f;
                    real tau_fCaSS = 80.0f / (1.0f + (actualCa_SS * actualCa_SS * 400.0f)) + 2.0f;

            #if defined(EPI) || defined(MCELL)
                    real s_inf = 1.0f / (1.0f + exp((actualV + 20.0f) / 5.0f));
                    real tau_s = 85.0f * exp(-(actualV + 45.0f) * (actualV + 45.0f) / 320.0f) + 5.0f / (1.0f + exp((actualV - 20.0f) / 5.0f)) + 3.0f;
            #endif // EPI || MCELL
            #ifdef ENDO
                    real s_inf = 1.0f / (1.0f + exp((actualV + 28.0f) / 5.0f));
                    real tau_s = 1000.0f * exp(-(actualV + 67.0f) * (actualV + 67.0f) / 1000.0f) + 8.0f;
            #endif // ENDO

                    real r_inf = 1.0f / (1.0f + exp((20.0f - actualV) / 6.0f));
                    real tau_r = 9.5f * exp(-(actualV + 40.0f) * (actualV + 40.0f) / 1800.0f) + 0.8f;

                    // Explicit method - auxiliary variables
                    real Ileak = V_leak * (actualCa_SR - actualCa_i);
                    real Iup = V_maxup / (1.0f + (K_up * K_up) / (actualCa_i * actualCa_i));
                    real k_CaSR = max_SR - (max_SR - min_SR) / (1.0f + (EC / actualCa_SR) * (EC / actualCa_SR));
                    real k1 = k1_prime / k_CaSR;
                    real k2 = k2_prime * k_CaSR;
                    real O = k1 * actualCa_SS * actualCa_SS * actualR_prime / (k3 + k1 * actualCa_SS * actualCa_SS);
                    real Irel = V_rel * O * (actualCa_SR - actualCa_SS);
                    real Ixfer = V_xfer * (actualCa_SS - actualCa_i);
                    real Ca_i_bufC; // !!!
                    Ca_i_bufC = 1.0f / (1.0f + bufC * K_bufC / ((K_bufC + actualCa_i) * (K_bufC + actualCa_i)));
                    real Ca_SR_bufSR; // !!!
                    Ca_SR_bufSR = 1.0f / (1.0f + bufSR * K_bufSR / ((K_bufSR + actualCa_SR) * (K_bufSR + actualCa_SR)));
                    real Ca_SS_bufSS; // !!!
                    Ca_SS_bufSS = 1.0f / (1.0f + bufSS * K_bufSS / ((K_bufSS + actualCa_SS) * (K_bufSS + actualCa_SS)));

                    // Explicit method - RHS of the state variables
                    real RHS_R_prime_term = (k4 * (1.0f - actualR_prime)) - k2 * actualCa_SS * actualR_prime;
                    real RHS_Ca_i_term = Ca_i_bufC * ((((Ileak - Iup) * V_SR / V_C) + Ixfer) - (((IbCa + IpCa) - 2.0f * INaCa) * Cm / (2.0f * V_C * F)));
                    real RHS_Ca_SR_term = Ca_SR_bufSR * (Iup - Ileak - Irel);
                    real RHS_Ca_SS_term = Ca_SS_bufSS * (((-ICaL * Cm / (2.0 * V_SS * F)) + (Irel * V_SR / V_SS)) - (Ixfer * V_C / V_SS));
                    real RHS_Na_i_term = -((INa + IbNa + 3.0f * INaK + 3.0f * INaCa) * Cm / (V_C * F));
                    real RHS_K_i_term = -((IK1 + Ito + IKr + IKs - 2.0f * INaK + IpK + stim) * Cm / (V_C * F));

                    // Rush-Larsen method - update approximations
                    real X_r1tilde = X_r1_inf - (X_r1_inf - actualX_r1) * exp(-(0.5f * delta_t) / (alpha_X_r1 * beta_X_r1));
                    real X_r2tilde = X_r2_inf - (X_r2_inf - actualX_r2) * exp(-(0.5f * delta_t) / (alpha_X_r2 * beta_X_r2));
                    real X_stilde = X_s_inf - (X_s_inf - actualX_s) * exp(-(0.5f * delta_t) / (alpha_X_s * beta_X_s + 80.0f));
                    real mtilde = m_inf - (m_inf - actualm) * exp(-(0.5f * delta_t) / (alpha_m * beta_m));
                    real htilde = h_inf - (h_inf - actualh) * exp(-(0.5f * delta_t) * (alpha_h + beta_h));
                    real jtilde = j_inf - (j_inf - actualj) * exp(-(0.5f * delta_t) * (alpha_j + beta_j));
                    real dtilde = inf - (inf - actuald) * exp(-(0.5f * delta_t) / (alpha_d * beta_d + gamma_d));
                    real ftilde = f_inf - (f_inf - actualf) * exp(-(0.5f * delta_t) / (alpha_f + beta_f + gamma_f));
                    real f2tilde = f2_inf - (f2_inf - actualf2) * exp(-(0.5f * delta_t) / (alpha_f2 + beta_f2 + gamma_f2));
                    real fCaSStilde = fCaSS_inf - (fCaSS_inf - actualfCaSS) * exp(-(0.5f * delta_t) / tau_fCaSS);
                    real stilde = s_inf - (s_inf - actuals) * exp(-(0.5f * delta_t) / tau_s);
                    real rtilde = r_inf - (r_inf - actualr) * exp(-(0.5f * delta_t) / tau_r);

                    // Explicit method - update approximations
                    real R_primetilde = actualR_prime + (0.5f * delta_t * RHS_R_prime_term);
                    real Ca_itilde = actualCa_i + (0.5f * delta_t * RHS_Ca_i_term);
                    real Ca_SRtilde = actualCa_SR + (0.5f * delta_t * RHS_Ca_SR_term);
                    real Ca_SStilde = actualCa_SS + (0.5f * delta_t * RHS_Ca_SS_term);
                    real Na_itilde = actualNa_i + (0.5f * delta_t * RHS_Na_i_term);
                    real K_itilde = actualK_i + (0.5f * delta_t * RHS_K_i_term);

                    // Auxiliary variables with actualVtilde
                    real VmENatilde = actualVtilde - (RTONF * log(Na_o / Na_itilde));
                    real E_Ktilde = (RTONF * log(K_o / K_itilde));
                    real VmEKtilde = actualVtilde - E_Ktilde;
                    real alpha_K1tilde = 0.1f / (1.0f + exp(0.06f * (actualVtilde - E_Ktilde - 200.0f)));
                    real beta_K1tilde = (3.0f * exp(0.0002f * (actualVtilde - E_Ktilde + 100.0f)) + exp(0.1f * (actualVtilde - E_Ktilde - 10.0f))) / (1.0f + exp(-0.5f * (actualVtilde - E_Ktilde)));
                    real E_Kstilde = (RTONF * log((K_o + p_KNa * Na_o) / (K_itilde + p_KNa * Na_itilde)));
                    real E_Catilde = 0.5f * RTONF * log(Ca_o / Ca_itilde);

                    // Currents with actualVtilde
                    real INatilde = G_Na * (mtilde * mtilde * mtilde) * htilde * jtilde * VmENatilde;
                    real IbNatilde = G_bNa * VmENatilde;
                    real IK1tilde = G_K1 * (alpha_K1tilde / (alpha_K1tilde + beta_K1tilde)) * VmEKtilde;
                    real Itotilde = G_to * rtilde * stilde * VmEKtilde;
                    real IKrtilde = G_Kr * sqrt(K_o / 5.4f) * X_r1tilde * X_r2tilde * VmEKtilde;
                    real IKstilde = G_Ks * X_stilde * X_stilde * (actualVtilde - E_Kstilde);
                    real ICaLtilde; // !!!
                    (actualVtilde < 15.0f - 1.0e-5f)
                        ? (ICaLtilde = G_CaL * dtilde * ftilde * f2tilde * fCaSStilde * 4.0f * (actualVtilde - 15.0f) * (F * F) * (0.25f * Ca_SStilde * exp(2.0f * (actualVtilde - 15.0f) * FONRT) - Ca_o) / (R * T * (exp(2.0f * (actualVtilde - 15.0f) * FONRT) - 1.0f)))
                        : (ICaLtilde = G_CaL * dtilde * ftilde * f2tilde * fCaSStilde * 2.0f * F * (0.25f * Ca_SStilde - Ca_o));
                    real INaKtilde = ((((p_KNa * K_o) / (K_o + K_mK)) * Na_itilde) / (Na_itilde + K_mNa)) / (1.0f + (0.1245f * exp(((-0.1f) * actualVtilde * FONRT))) + (0.0353f * exp(((-actualVtilde) * FONRT))));
                    real INaCatilde; // !!!
                    INaCatilde = (k_NaCa * ((exp((gamma_I_NaCa * actualVtilde * FONRT)) * (Na_itilde * Na_itilde * Na_itilde) * Ca_o) - (exp(((gamma_I_NaCa - 1.0f) * actualVtilde * FONRT)) * (Na_o * Na_o * Na_o) * Ca_itilde * alpha))) / (((K_mNa_i * K_mNa_i * K_mNa_i) + (Na_o * Na_o * Na_o)) * (K_mCa + Ca_o) * (1.0f + (k_sat * exp(((gamma_I_NaCa)*actualVtilde * FONRT)))));
                    real IpCatilde = (G_pCa * Ca_itilde) / (K_pCa + Ca_itilde);
                    real IpKtilde = (G_pK * VmEKtilde) / (1.0f + exp((25.0f - actualVtilde) / 5.98f));
                    real IbCatilde = G_bCa * (actualVtilde - E_Catilde);

                    // RHS of the main equation with actualVtilde
                    real RHS_Vtilde_term = INatilde + IbNatilde + IK1tilde + Itotilde + IKrtilde + IKstilde + ICaLtilde + INaKtilde + INaCatilde + IpCatilde + IpKtilde + IbCatilde;
                    partRHS[i][j] = delta_t * (stim - RHS_Vtilde_term);

                    // Update state variables
                    // RHS of the state variables with tilde approximations
                    // Rush-Larsen method - auxiliary variables
                    real X_r1_inftilde = 1.0f / (1.0f + exp((-26.0f - actualVtilde) / 7.0f));
                    real alpha_X_r1tilde = 450.0f / (1.0f + exp((-45.0f - actualVtilde) / 10.0f));
                    real beta_X_r1tilde = 6.0f / (1.0f + exp((30.0f + actualVtilde) / 11.5f));

                    real X_r2_inftilde = 1.0f / (1.0f + exp((actualVtilde + 88.0f) / 24.0f));
                    real alpha_X_r2tilde = 3.0f / (1.0f + exp((-60.0f - actualVtilde) / 20.0f));
                    real beta_X_r2tilde = 1.12f / (1.0f + exp((actualVtilde - 60.0f) / 20.0f));

                    real X_s_inftilde = 1.0f / (1.0f + exp((-5.0f - actualVtilde) / 14.0f));
                    real alpha_X_stilde = 1400.0f / sqrt(1.0f + exp((5.0f - actualVtilde) / 6.0f));
                    real beta_X_stilde = 1.0f / (1.0f + exp((-35.0f + actualVtilde) / 15.0f));

                    real m_inftilde = 1.0f / ((1.0f + exp((-56.86 - actualVtilde) / 9.03f)) * (1.0f + exp((-56.86f - actualVtilde) / 9.03f)));
                    real alpha_mtilde = 1.0f / (1.0f + exp((-60.0f - actualVtilde) / 5.0f));
                    real beta_mtilde = 0.1f / (1.0f + exp((actualVtilde + 35.0f) / 5.0f)) + (0.1f / (1.0f + exp((actualVtilde - 50.0f) / 200.0f)));

                    real h_inftilde = 1.0f / ((1.0f + exp((actualVtilde + 71.55f) / 7.43f)) * (1.0f + exp((actualVtilde + 71.55f) / 7.43f)));
                    real alpha_htilde;
                    (actualVtilde < -40.0f)
                        ? (alpha_htilde = 0.057f * exp(-(actualVtilde + 80.0f) / 6.8f))
                        : (alpha_htilde = 0.0f);
                    real beta_htilde;
                    (actualVtilde < -40.0f)
                        ? (beta_htilde = 2.7f * exp(0.079f * actualVtilde) + 3.1f * 1.0e5f * exp(0.3485f * actualVtilde))
                        : (beta_htilde = 0.77f / (0.13f * (1.0f + exp((actualVtilde + 10.66f) / -11.1f))));

                    real j_inftilde = 1.0f / ((1.0f + exp((actualVtilde + 71.55f) / 7.43f)) * (1.0f + exp((actualVtilde + 71.55f) / 7.43f)));
                    real alpha_jtilde;
                    (actualVtilde < -40.0f)
                        ? (alpha_jtilde = ((-25428.0f * exp(0.2444f * actualVtilde) - (6.948e-6f * exp((-0.04391f) * actualVtilde))) * (actualVtilde + 37.78f)) / (1.0f + exp(0.311f * (actualVtilde + 79.23f))))
                        : (alpha_jtilde = 0.0f);
                    real beta_jtilde;
                    (actualVtilde < -40.0f)
                        ? (beta_jtilde = (0.02424f * exp(-0.01052f * actualVtilde)) / (1.0f + exp(-0.1378f * (actualVtilde + 40.14f))))
                        : (beta_jtilde = (0.6f * exp(0.057f * actualVtilde)) / (1.0f + exp(-0.1f * (actualVtilde + 32.0f))));

                    real d_inftilde = 1.0f / (1.0f + exp((-8.0f - actualVtilde) / 7.5f));
                    real alpha_dtilde = 1.4f / (1.0f + exp((-35.0f - actualVtilde) / 13.0f)) + 0.25f;
                    real beta_dtilde = 1.4f / (1.0f + exp((actualVtilde + 5.0f) / 5.0f));
                    real gamma_dtilde = 1.0f / (1.0f + exp((50.0f - actualVtilde) / 20.0f));

                    real f_inftilde = 1.0f / (1.0f + exp((actualVtilde + 20.0f) / 7.0f));
                    real alpha_ftilde = 1102.5f * exp(-(actualVtilde + 27.0f) * (actualVtilde + 27.0f) / 225.0f);
                    real beta_ftilde = 200.0f / (1.0f + exp((13.0f - actualVtilde) / 10.0f));
                    real gamma_ftilde = 180.0f / (1.0f + exp((actualVtilde + 30.0f) / 10.0f)) + 20.0f;

                    real f2_inftilde = 0.67f / (1.0f + exp((actualVtilde + 35.0f) / 7.0f)) + 0.33f;
                    real alpha_f2tilde; // !!!
                    alpha_f2tilde = 562.0f * exp(-(actualVtilde + 27.0f) * (actualVtilde + 27.0f) / 240.0f);
                    real beta_f2tilde = 31.0f / (1.0f + exp((25.0f - actualVtilde) / 10.0f));
                    real gamma_f2tilde; // !!!
                    gamma_f2tilde = 80.0f / (1.0f + exp((30.0f + actualVtilde) / 10.0f));

                    real fCaSS_inftilde = 0.6f / (1.0f + (Ca_SStilde * Ca_SStilde * 400.0f)) + 0.4f;
                    real tau_fCaSStilde = 80.0f / (1.0f + (Ca_SStilde * Ca_SStilde * 400.0f)) + 2.0f;

            #if defined(EPI) || defined(MCELL)
                    real s_inftilde = 1.0f / (1.0f + exp((actualVtilde + 20.0f) / 5.0f));
                    real tau_stilde = 85.0f * exp(-(actualVtilde + 45.0f) * (actualVtilde + 45.0f) / 320.0f) + 5.0f / (1.0f + exp((actualVtilde - 20.0f) / 5.0f)) + 3.0f;
            #endif // EPI || MCELL
            #ifdef ENDO
                    real s_inftilde = 1.0f / (1.0f + exp((actualVtilde + 28.0f) / 5.0f));
                    real tau_stilde = 1000.0f * exp(-(actualVtilde + 67.0f) * (actualVtilde + 67.0f) / 1000.0f) + 8.0f;
            #endif // ENDO

                    real r_inftilde = 1.0f / (1.0f + exp((20.0f - actualVtilde) / 6.0f));
                    real tau_rtilde = 9.5f * exp(-(actualVtilde + 40.0f) * (actualVtilde + 40.0f) / 1800.0f) + 0.8f;

                    // Explicit method - auxiliary variables
                    real Ileaktilde = V_leak * (Ca_SRtilde - Ca_itilde);
                    real Iuptilde = V_maxup / (1.0f + (K_up * K_up) / (Ca_itilde * Ca_itilde));
                    real k_CaSRtilde = max_SR - (max_SR - min_SR) / (1.0f + (EC / Ca_SRtilde) * (EC / Ca_SRtilde));
                    real k1tilde = k1_prime / k_CaSRtilde;
                    real k2tilde = k2_prime * k_CaSRtilde;
                    real Otilde = k1tilde * Ca_SStilde * Ca_SStilde * R_primetilde / (k3 + k1tilde * Ca_SStilde * Ca_SStilde);
                    real Ireltilde = V_rel * Otilde * (Ca_SRtilde - Ca_SStilde);
                    real Ixfertilde = V_xfer * (Ca_SStilde - Ca_itilde);
                    real Ca_i_bufCtilde; // !!!
                    Ca_i_bufCtilde = 1.0f / (1.0f + bufC * K_bufC / ((K_bufC + Ca_itilde) * (K_bufC + Ca_itilde)));
                    real Ca_SR_bufSRtilde; // !!!
                    Ca_SR_bufSRtilde = 1.0f / (1.0f + bufSR * K_bufSR / ((K_bufSR + Ca_SRtilde) * (K_bufSR + Ca_SRtilde)));
                    real Ca_SS_bufSStilde; // !!!
                    Ca_SS_bufSStilde = 1.0f / (1.0f + bufSS * K_bufSS / ((K_bufSS + Ca_SStilde) * (K_bufSS + Ca_SStilde)));

                    // Explicit method - RHS of the state variables
                    real RHS_R_primetilde_term = (k4 * (1.0f - R_primetilde)) - k2tilde * Ca_SStilde * R_primetilde;
                    real RHS_Ca_itilde_term = Ca_i_bufCtilde * ((((Ileaktilde - Iuptilde) * V_SR / V_C) + Ixfertilde) - (((IbCatilde + IpCatilde) - 2.0f * INaCatilde) * Cm / (2.0f * V_C * F)));
                    real RHS_Ca_SRtilde_term = Ca_SR_bufSRtilde * (Iuptilde - Ileaktilde - Ireltilde);
                    real RHS_Ca_SStilde_term = Ca_SS_bufSStilde * (((-ICaLtilde * Cm / (2.0 * V_SS * F)) + (Ireltilde * V_SR / V_SS)) - (Ixfertilde * V_C / V_SS));
                    real RHS_Na_itilde_term = -((INatilde + IbNatilde + 3.0f * INaKtilde + 3.0f * INaCatilde) * Cm / (V_C * F));
                    real RHS_K_itilde_term = -((IK1tilde + Itotilde + IKrtilde + IKstilde - 2.0f * INaKtilde + IpKtilde + stim) * Cm / (V_C * F));

                    // Rush-Larsen method - update approximations
                    X_r1[i][j] = X_r1_inftilde - (X_r1_inftilde - actualX_r1) * exp(-delta_t / (alpha_X_r1tilde * beta_X_r1tilde));
                    X_r2[i][j] = X_r2_inftilde - (X_r2_inftilde - actualX_r2) * exp(-delta_t / (alpha_X_r2tilde * beta_X_r2tilde));
                    X_s[i][j] = X_s_inftilde - (X_s_inftilde - actualX_s) * exp(-delta_t / (alpha_X_stilde * beta_X_stilde + 80.0f));
                    m[i][j] = m_inftilde - (m_inftilde - actualm) * exp(-delta_t / (alpha_mtilde * beta_mtilde));
                    h[i][j] = h_inftilde - (h_inftilde - actualh) * exp(-delta_t * (alpha_htilde + beta_htilde));
                    _j[i][j] = j_inftilde - (j_inftilde - actualj) * exp(-delta_t * (alpha_jtilde + beta_jtilde));
                    d[i][j] = d_inftilde - (d_inftilde - actuald) * exp(-delta_t / (alpha_dtilde * beta_dtilde + gamma_dtilde));
                    f[i][j] = f_inftilde - (f_inftilde - actualf) * exp(-delta_t / (alpha_ftilde + beta_ftilde + gamma_ftilde));
                    f2[i][j] = f2_inftilde - (f2_inftilde - actualf2) * exp(-delta_t / (alpha_f2tilde + beta_f2tilde + gamma_f2tilde));
                    fCaSS[i][j] = fCaSS_inftilde - (fCaSS_inftilde - actualfCaSS) * exp(-delta_t / tau_fCaSStilde);
                    s[i][j] = s_inftilde - (s_inftilde - actuals) * exp(-delta_t / tau_stilde);
                    r[i][j] = r_inftilde - (r_inftilde - actualr) * exp(-delta_t / tau_rtilde);

                    // Explicit method - update approximations
                    R_prime[i][j] = actualR_prime + (delta_t * RHS_R_primetilde_term);
                    Ca_i[i][j] = actualCa_i + (delta_t * RHS_Ca_itilde_term);
                    Ca_SR[i][j] = actualCa_SR + (delta_t * RHS_Ca_SRtilde_term);
                    Ca_SS[i][j] = actualCa_SS + (delta_t * RHS_Ca_SStilde_term);
                    Na_i[i][j] = actualNa_i + (delta_t * RHS_Na_itilde_term);
                    K_i[i][j] = actualK_i + (delta_t * RHS_K_itilde_term);
#endif // TT2
#endif // MONODOMAIN
                }
            }

            // ================================================!
            //  Calculate V at n+1/2 -> Result goes to RHS     !
            // ================================================!
            for (int j = 0; j < N; j++)
            {
                for (int i = 0; i < N; i++)
                {
                    real actualV = V[i][j];
                    real tau = 0.5f;
                    if (strcmp(method, "theta-ADI") == 0)
                        tau =  1.0f - theta;
                    
                    diff_term =  diff_coeff * tau * phi * (V[i][lim(j - 1, N)] - 2.0f * actualV + V[i][lim(j + 1, N)]);
                    LS_b[i] = actualV + diff_term + 0.5f * partRHS[i][j];
                }

                tridiag(la, lb, lc, c_prime, d_prime, N, LS_b, result);
                for (int i = 0; i < N; i++)
                {
                    RHS[i][j] = result[i];
                }
            }

            // ================================================!
            //  Calculate V at n+1 -> Result goes to V         !
            // ================================================!
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    real actualV = RHS[i][j];
                    real tau = 0.5f;
                    if (strcmp(method, "theta-ADI") == 0)
                        tau =  1.0f - theta;
                    
                    diff_term =  diff_coeff * tau * phi * (RHS[lim(i - 1, N)][j] - 2.0f * actualV + RHS[lim(i + 1, N)][j]);
                    LS_b[j] = actualV + diff_term + 0.5f * partRHS[i][j];
                }

                tridiag(la, lb, lc, c_prime, d_prime, N, LS_b, result);
                for (int j = 0; j < N; j++)
                {
                    V[i][j] = result[j];
                }
            }

#ifdef SAVE_FRAMES
            // If save frames is true and time step is multiple of frame save rate
            if (timeStepCounter % frameSaveRate == 0)
            {
                // Save frame
                saveFrame(fpFrames, actualTime, V, N);
                printf("Frame at time %.2lf ms saved to %s\n", actualTime, framesPath);
            }
#endif // SAVE_FRAMES
#ifndef CONVERGENCE_ANALYSIS
            // Calculate stim velocity
            if (!stim_velocity_measured)
            {
                real begin = L / 3.0f;
                real end = 2.0f * begin;

                if (!aux_stim_velocity_flag)
                {
                    int first_point_index = round(begin / delta_x) + 1;
                    if (V[0][first_point_index] > 10.0f)
                    {
                        first_point_time = actualTime;
                        aux_stim_velocity_flag = true;
                    }
                }
                else
                {
                    int last_point_index = round(end / delta_x) + 1;
                    if (V[0][last_point_index] > 10.0f)
                    {
                        last_point_time = actualTime;
                        stim_velocity = (end - begin) / (last_point_time - first_point_time); // cm/ms
                        stim_velocity = stim_velocity * 10.0f; // m/s
                        stim_velocity_measured = true;
                        printf("Stim velocity (measured from %.2f to %.2f cm) is %lf m/s\n", begin, end, stim_velocity);
                    }
                }
            }
#endif // not CONVERGENCE_ANALYSIS

            // Update time step counter
            timeStepCounter++;
        }
    }
#endif // LINMONO || MONODOMAIN

#else // if not def CABLEEQ
    if (strcmp(method, "theta-RK2") == 0)
    {
        while (timeStepCounter < M)
        {
            // Get time step
            actualTime = time[timeStepCounter];

            #ifdef CABLEEQ
            // Get info for Action Potential
            AP[timeStepCounter] = V[APCellIndex]; 
            #endif // CABLEEQ

            // ================================================!
            //  Calcula Approx.                                !
            // ================================================!
            real diff_term = 0.0f;
            for (int i = 0; i < N; i++)
            {
#ifdef AFHN
                real actualV = V[i];
                real actualW = W[i];

                diff_term = diff_coeff * phi * (V[lim(i - 1, N)] - 2.0f * actualV + V[lim(i + 1, N)]);
                real RHS_V_term = RHS_V(actualV, actualW) / (Cm * chi);

                // Stimulation
                real stim = 0.0f;
                for (int si = 0; si < numberOfStimuli; si++)
                {
                    if (actualTime >= stimuli[si].begin && actualTime <= stimuli[si].begin + stimuli[si].duration && i >= stimuli[si].xMinDisc && i <= stimuli[si].xMaxDisc)
                    {
                        stim = stimuli[si].strength;
                        break;
                    }
                }

                // Calculate Vtilde -> utilde = u^n + 0.5 * dt * (A*u^n + R(u^n))
                real actualVtilde = actualV + 0.5f * diff_term + (0.5f * delta_t * (stim - RHS_V_term));

                // Calculate W approximation
                real RHS_W_term = RHS_W(actualV, actualW);
                real Wtilde = actualW + (0.5f * delta_t * RHS_W_term);

                // Preparing part of the RHS of the following linear systems
                real RHS_Vtilde_term = RHS_V(actualVtilde, Wtilde) / (Cm * chi);
                partRHS[i] = delta_t * (stim - RHS_Vtilde_term);

                // Update Wn+1 with RK2 -> Wn+1 = Wn + dt*R(V*, W*)
                W[i] = actualW + delta_t * RHS_W(actualVtilde, Wtilde);
#endif // AFHN
#ifdef TT2
                real actualV = V[i];
                real im1V = V[lim(i - 1, N)];
                real ip1V = V[lim(i + 1, N)];

                diff_term = diff_coeff * phi * (im1V - 2.0f * actualV + ip1V);

                // State variables
                real actualX_r1 = X_r1[i];
                real actualX_r2 = X_r2[i];
                real actualX_s = X_s[i];
                real actualm = m[i];
                real actualh = h[i];
                real actualj = _j[i];
                real actuald = d[i];
                real actualf = f[i];
                real actualf2 = f2[i];
                real actualfCaSS = fCaSS[i];
                real actuals = s[i];
                real actualr = r[i];
                real actualCa_i = Ca_i[i];
                real actualCa_SR = Ca_SR[i];
                real actualCa_SS = Ca_SS[i];
                real actualR_prime = R_prime[i];
                real actualNa_i = Na_i[i];
                real actualK_i = K_i[i];

                // Auxiliary variables
                real VmENa = actualV - (RTONF * log(Na_o / actualNa_i));
                real E_K = (RTONF * log(K_o / actualK_i));
                real VmEK = actualV - E_K;
                real alpha_K1 = 0.1f / (1.0f + exp(0.06f * (actualV - E_K - 200.0f)));
                real beta_K1 = (3.0f * exp(0.0002f * (actualV - E_K + 100.0f)) + exp(0.1f * (actualV - E_K - 10.0f))) / (1.0f + exp(-0.5f * (actualV - E_K)));
                real E_Ks = (RTONF * log((K_o + p_KNa * Na_o) / (actualK_i + p_KNa * actualNa_i)));
                real E_Ca = 0.5f * RTONF * log(Ca_o / actualCa_i);

                // Currents
                real INa = G_Na * (actualm * actualm * actualm) * actualh * actualj * VmENa;
                real IbNa = G_bNa * VmENa;
                real IK1 = G_K1 * (alpha_K1 / (alpha_K1 + beta_K1)) * VmEK;
                real Ito = G_to * actualr * actuals * VmEK;
                real IKr = G_Kr * sqrt(K_o / 5.4f) * actualX_r1 * actualX_r2 * VmEK;
                real IKs = G_Ks * actualX_s * actualX_s * (actualV - E_Ks);
                real ICaL; // !!!
                (actualV < 15.0f - 1.0e-5f)
                    ? (ICaL = G_CaL * actuald * actualf * actualf2 * actualfCaSS * 4.0f * (actualV - 15.0f) * (F * F) * (0.25f * actualCa_SS * exp(2.0f * (actualV - 15.0f) * FONRT) - Ca_o) / (R * T * (exp(2.0f * (actualV - 15.0f) * FONRT) - 1.0f)))
                    : (ICaL = G_CaL * actuald * actualf * actualf2 * actualfCaSS * 2.0f * F * (0.25f * actualCa_SS - Ca_o));
                real INaK = ((((p_KNa * K_o) / (K_o + K_mK)) * actualNa_i) / (actualNa_i + K_mNa)) / (1.0f + (0.1245f * exp(((-0.1f) * actualV * FONRT))) + (0.0353f * exp(((-actualV) * FONRT))));
                real INaCa; // !!!
                INaCa = (k_NaCa * ((exp((gamma_I_NaCa * actualV * FONRT)) * (actualNa_i * actualNa_i * actualNa_i) * Ca_o) - (exp(((gamma_I_NaCa - 1.0f) * actualV * FONRT)) * (Na_o * Na_o * Na_o) * actualCa_i * alpha))) / (((K_mNa_i * K_mNa_i * K_mNa_i) + (Na_o * Na_o * Na_o)) * (K_mCa + Ca_o) * (1.0f + (k_sat * exp(((gamma_I_NaCa)*actualV * FONRT)))));
                real IpCa = (G_pCa * actualCa_i) / (K_pCa + actualCa_i);
                real IpK = (G_pK * VmEK) / (1.0f + exp((25.0f - actualV) / 5.98f));
                real IbCa = G_bCa * (actualV - E_Ca);

                // RHS of the main equation
                real RHS_V_term = INa + IbNa + IK1 + Ito + IKr + IKs + ICaL + INaK + INaCa + IpCa + IpK + IbCa;

                real stim = 0.0f;
                for (int si = 0; si < numberOfStimuli; si++)
                {
                    if (actualTime >= stimuli[si].begin && actualTime <= stimuli[si].begin + stimuli[si].duration && i >= stimuli[si].xMinDisc && i <= stimuli[si].xMaxDisc)
                    {
                        stim = stimuli[si].strength;
                        break;
                    }
                }

                // Calculate Vtilde -> utilde = u^n + 0.5 * dt * (A*u^n + R(u^n))
                real actualVtilde = actualV + 0.5f * diff_term + (0.5f * delta_t * (stim - RHS_V_term));

                // Preparing part of the RHS of the following linear systems
                // Calculate approximation for state variables
                // Rush-Larsen method - auxiliary variables
                real X_r1_inf = 1.0f / (1.0f + exp((-26.0f - actualV) / 7.0f));
                real alpha_X_r1 = 450.0f / (1.0f + exp((-45.0f - actualV) / 10.0f));
                real beta_X_r1 = 6.0f / (1.0f + exp((30.0f + actualV) / 11.5f));

                real X_r2_inf = 1.0f / (1.0f + exp((actualV + 88.0f) / 24.0f));
                real alpha_X_r2 = 3.0f / (1.0f + exp((-60.0f - actualV) / 20.0f));
                real beta_X_r2 = 1.12f / (1.0f + exp((actualV - 60.0f) / 20.0f));

                real X_s_inf = 1.0f / (1.0f + exp((-5.0f - actualV) / 14.0f));
                real alpha_X_s = 1400.0f / sqrt(1.0f + exp((5.0f - actualV) / 6.0f));
                real beta_X_s = 1.0f / (1.0f + exp((-35.0f + actualV) / 15.0f));

                real m_inf = 1.0f / ((1.0f + exp((-56.86 - actualV) / 9.03f)) * (1.0f + exp((-56.86f - actualV) / 9.03f)));
                real alpha_m = 1.0f / (1.0f + exp((-60.0f - actualV) / 5.0f));
                real beta_m = 0.1f / (1.0f + exp((actualV + 35.0f) / 5.0f)) + (0.1f / (1.0f + exp((actualV - 50.0f) / 200.0f)));

                real h_inf = 1.0f / ((1.0f + exp((actualV + 71.55f) / 7.43f)) * (1.0f + exp((actualV + 71.55f) / 7.43f)));
                real alpha_h;
                (actualV < -40.0f)
                    ? (alpha_h = 0.057f * exp(-(actualV + 80.0f) / 6.8f))
                    : (alpha_h = 0.0f);
                real beta_h;
                (actualV < -40.0f)
                    ? (beta_h = 2.7f * exp(0.079f * actualV) + 3.1f * 1.0e5f * exp(0.3485f * actualV))
                    : (beta_h = 0.77f / (0.13f * (1.0f + exp((actualV + 10.66f) / -11.1f))));

                real j_inf = 1.0f / ((1.0f + exp((actualV + 71.55f) / 7.43f)) * (1.0f + exp((actualV + 71.55f) / 7.43f)));
                real alpha_j;
                (actualV < -40.0f)
                    ? (alpha_j = ((-25428.0f * exp(0.2444f * actualV) - (6.948e-6f * exp((-0.04391f) * actualV))) * (actualV + 37.78f)) / (1.0f + exp(0.311f * (actualV + 79.23f))))
                    : (alpha_j = 0.0f);
                real beta_j;
                (actualV < -40.0f)
                    ? (beta_j = (0.02424f * exp(-0.01052f * actualV)) / (1.0f + exp(-0.1378f * (actualV + 40.14f))))
                    : (beta_j = (0.6f * exp(0.057f * actualV)) / (1.0f + exp(-0.1f * (actualV + 32.0f))));

                real inf = 1.0f / (1.0f + exp((-8.0f - actualV) / 7.5f));
                real alpha_d = 1.4f / (1.0f + exp((-35.0f - actualV) / 13.0f)) + 0.25f;
                real beta_d = 1.4f / (1.0f + exp((actualV + 5.0f) / 5.0f));
                real gamma_d = 1.0f / (1.0f + exp((50.0f - actualV) / 20.0f));

                real f_inf = 1.0f / (1.0f + exp((actualV + 20.0f) / 7.0f));
                real alpha_f = 1102.5f * exp(-(actualV + 27.0f) * (actualV + 27.0f) / 225.0f);
                real beta_f = 200.0f / (1.0f + exp((13.0f - actualV) / 10.0f));
                real gamma_f = 180.0f / (1.0f + exp((actualV + 30.0f) / 10.0f)) + 20.0f;

                real f2_inf = 0.67f / (1.0f + exp((actualV + 35.0f) / 7.0f)) + 0.33f;
                real alpha_f2; // !!!
                alpha_f2 = 562.0f * exp(-(actualV + 27.0f) * (actualV + 27.0f) / 240.0f);
                real beta_f2 = 31.0f / (1.0f + exp((25.0f - actualV) / 10.0f));
                real gamma_f2; // !!!
                gamma_f2 = 80.0f / (1.0f + exp((30.0f + actualV) / 10.0f));

                real fCaSS_inf = 0.6f / (1.0f + (actualCa_SS * actualCa_SS * 400.0f)) + 0.4f;
                real tau_fCaSS = 80.0f / (1.0f + (actualCa_SS * actualCa_SS * 400.0f)) + 2.0f;

        #if defined(EPI) || defined(MCELL)
                real s_inf = 1.0f / (1.0f + exp((actualV + 20.0f) / 5.0f));
                real tau_s = 85.0f * exp(-(actualV + 45.0f) * (actualV + 45.0f) / 320.0f) + 5.0f / (1.0f + exp((actualV - 20.0f) / 5.0f)) + 3.0f;
        #endif // EPI || MCELL
        #ifdef ENDO
                real s_inf = 1.0f / (1.0f + exp((actualV + 28.0f) / 5.0f));
                real tau_s = 1000.0f * exp(-(actualV + 67.0f) * (actualV + 67.0f) / 1000.0f) + 8.0f;
        #endif // ENDO

                real r_inf = 1.0f / (1.0f + exp((20.0f - actualV) / 6.0f));
                real tau_r = 9.5f * exp(-(actualV + 40.0f) * (actualV + 40.0f) / 1800.0f) + 0.8f;

                // Explicit method - auxiliary variables
                real Ileak = V_leak * (actualCa_SR - actualCa_i);
                real Iup = V_maxup / (1.0f + (K_up * K_up) / (actualCa_i * actualCa_i));
                real k_CaSR = max_SR - (max_SR - min_SR) / (1.0f + (EC / actualCa_SR) * (EC / actualCa_SR));
                real k1 = k1_prime / k_CaSR;
                real k2 = k2_prime * k_CaSR;
                real O = k1 * actualCa_SS * actualCa_SS * actualR_prime / (k3 + k1 * actualCa_SS * actualCa_SS);
                real Irel = V_rel * O * (actualCa_SR - actualCa_SS);
                real Ixfer = V_xfer * (actualCa_SS - actualCa_i);
                real Ca_i_bufC; // !!!
                Ca_i_bufC = 1.0f / (1.0f + bufC * K_bufC / ((K_bufC + actualCa_i) * (K_bufC + actualCa_i)));
                real Ca_SR_bufSR; // !!!
                Ca_SR_bufSR = 1.0f / (1.0f + bufSR * K_bufSR / ((K_bufSR + actualCa_SR) * (K_bufSR + actualCa_SR)));
                real Ca_SS_bufSS; // !!!
                Ca_SS_bufSS = 1.0f / (1.0f + bufSS * K_bufSS / ((K_bufSS + actualCa_SS) * (K_bufSS + actualCa_SS)));

                // Explicit method - RHS of the state variables
                real RHS_R_prime_term = (k4 * (1.0f - actualR_prime)) - k2 * actualCa_SS * actualR_prime;
                real RHS_Ca_i_term = Ca_i_bufC * ((((Ileak - Iup) * V_SR / V_C) + Ixfer) - (((IbCa + IpCa) - 2.0f * INaCa) * Cm / (2.0f * V_C * F)));
                real RHS_Ca_SR_term = Ca_SR_bufSR * (Iup - Ileak - Irel);
                real RHS_Ca_SS_term = Ca_SS_bufSS * (((-ICaL * Cm / (2.0 * V_SS * F)) + (Irel * V_SR / V_SS)) - (Ixfer * V_C / V_SS));
                real RHS_Na_i_term = -((INa + IbNa + 3.0f * INaK + 3.0f * INaCa) * Cm / (V_C * F));
                real RHS_K_i_term = -((IK1 + Ito + IKr + IKs - 2.0f * INaK + IpK + stim) * Cm / (V_C * F));

                // Rush-Larsen method - update approximations
                real X_r1tilde = X_r1_inf - (X_r1_inf - actualX_r1) * exp(-(0.5f * delta_t) / (alpha_X_r1 * beta_X_r1));
                real X_r2tilde = X_r2_inf - (X_r2_inf - actualX_r2) * exp(-(0.5f * delta_t) / (alpha_X_r2 * beta_X_r2));
                real X_stilde = X_s_inf - (X_s_inf - actualX_s) * exp(-(0.5f * delta_t) / (alpha_X_s * beta_X_s + 80.0f));
                real mtilde = m_inf - (m_inf - actualm) * exp(-(0.5f * delta_t) / (alpha_m * beta_m));
                real htilde = h_inf - (h_inf - actualh) * exp(-(0.5f * delta_t) * (alpha_h + beta_h));
                real jtilde = j_inf - (j_inf - actualj) * exp(-(0.5f * delta_t) * (alpha_j + beta_j));
                real dtilde = inf - (inf - actuald) * exp(-(0.5f * delta_t) / (alpha_d * beta_d + gamma_d));
                real ftilde = f_inf - (f_inf - actualf) * exp(-(0.5f * delta_t) / (alpha_f + beta_f + gamma_f));
                real f2tilde = f2_inf - (f2_inf - actualf2) * exp(-(0.5f * delta_t) / (alpha_f2 + beta_f2 + gamma_f2));
                real fCaSStilde = fCaSS_inf - (fCaSS_inf - actualfCaSS) * exp(-(0.5f * delta_t) / tau_fCaSS);
                real stilde = s_inf - (s_inf - actuals) * exp(-(0.5f * delta_t) / tau_s);
                real rtilde = r_inf - (r_inf - actualr) * exp(-(0.5f * delta_t) / tau_r);

                // Explicit method - update approximations
                real R_primetilde = actualR_prime + (0.5f * delta_t * RHS_R_prime_term);
                real Ca_itilde = actualCa_i + (0.5f * delta_t * RHS_Ca_i_term);
                real Ca_SRtilde = actualCa_SR + (0.5f * delta_t * RHS_Ca_SR_term);
                real Ca_SStilde = actualCa_SS + (0.5f * delta_t * RHS_Ca_SS_term);
                real Na_itilde = actualNa_i + (0.5f * delta_t * RHS_Na_i_term);
                real K_itilde = actualK_i + (0.5f * delta_t * RHS_K_i_term);

                // Auxiliary variables with Vtilde
                real VmENatilde = actualVtilde - (RTONF * log(Na_o / Na_itilde));
                real E_Ktilde = (RTONF * log(K_o / K_itilde));
                real VmEKtilde = actualVtilde - E_Ktilde;
                real alpha_K1tilde = 0.1f / (1.0f + exp(0.06f * (actualVtilde - E_Ktilde - 200.0f)));
                real beta_K1tilde = (3.0f * exp(0.0002f * (actualVtilde - E_Ktilde + 100.0f)) + exp(0.1f * (actualVtilde - E_Ktilde - 10.0f))) / (1.0f + exp(-0.5f * (actualVtilde - E_Ktilde)));
                real E_Kstilde = (RTONF * log((K_o + p_KNa * Na_o) / (K_itilde + p_KNa * Na_itilde)));
                real E_Catilde = 0.5f * RTONF * log(Ca_o / Ca_itilde);

                // Currents with Vtilde
                real INatilde = G_Na * (mtilde * mtilde * mtilde) * htilde * jtilde * VmENatilde;
                real IbNatilde = G_bNa * VmENatilde;
                real IK1tilde = G_K1 * (alpha_K1tilde / (alpha_K1tilde + beta_K1tilde)) * VmEKtilde;
                real Itotilde = G_to * rtilde * stilde * VmEKtilde;
                real IKrtilde = G_Kr * sqrt(K_o / 5.4f) * X_r1tilde * X_r2tilde * VmEKtilde;
                real IKstilde = G_Ks * X_stilde * X_stilde * (actualVtilde - E_Kstilde);
                real ICaLtilde; // !!!
                (actualVtilde < 15.0f - 1.0e-5f)
                    ? (ICaLtilde = G_CaL * dtilde * ftilde * f2tilde * fCaSStilde * 4.0f * (actualVtilde - 15.0f) * (F * F) * (0.25f * Ca_SStilde * exp(2.0f * (actualVtilde - 15.0f) * FONRT) - Ca_o) / (R * T * (exp(2.0f * (actualVtilde - 15.0f) * FONRT) - 1.0f)))
                    : (ICaLtilde = G_CaL * dtilde * ftilde * f2tilde * fCaSStilde * 2.0f * F * (0.25f * Ca_SStilde - Ca_o));
                real INaKtilde = ((((p_KNa * K_o) / (K_o + K_mK)) * Na_itilde) / (Na_itilde + K_mNa)) / (1.0f + (0.1245f * exp(((-0.1f) * actualVtilde * FONRT))) + (0.0353f * exp(((-actualVtilde) * FONRT))));
                real INaCatilde; // !!!
                INaCatilde = (k_NaCa * ((exp((gamma_I_NaCa * actualVtilde * FONRT)) * (Na_itilde * Na_itilde * Na_itilde) * Ca_o) - (exp(((gamma_I_NaCa - 1.0f) * actualVtilde * FONRT)) * (Na_o * Na_o * Na_o) * Ca_itilde * alpha))) / (((K_mNa_i * K_mNa_i * K_mNa_i) + (Na_o * Na_o * Na_o)) * (K_mCa + Ca_o) * (1.0f + (k_sat * exp(((gamma_I_NaCa)*actualVtilde * FONRT)))));
                real IpCatilde = (G_pCa * Ca_itilde) / (K_pCa + Ca_itilde);
                real IpKtilde = (G_pK * VmEKtilde) / (1.0f + exp((25.0f - actualVtilde) / 5.98f));
                real IbCatilde = G_bCa * (actualVtilde - E_Catilde);

                // part of RHS of the main equation with Vtilde
                real RHS_Vtilde_term = INatilde + IbNatilde + IK1tilde + Itotilde + IKrtilde + IKstilde + ICaLtilde + INaKtilde + INaCatilde + IpCatilde + IpKtilde + IbCatilde;
                partRHS[i] = delta_t * (stim - RHS_Vtilde_term);

                // Update state variables
                // RHS of the state variables with tilde approximations
                // Rush-Larsen method - auxiliary variables
                real X_r1_inftilde = 1.0f / (1.0f + exp((-26.0f - actualVtilde) / 7.0f));
                real alpha_X_r1tilde = 450.0f / (1.0f + exp((-45.0f - actualVtilde) / 10.0f));
                real beta_X_r1tilde = 6.0f / (1.0f + exp((30.0f + actualVtilde) / 11.5f));

                real X_r2_inftilde = 1.0f / (1.0f + exp((actualVtilde + 88.0f) / 24.0f));
                real alpha_X_r2tilde = 3.0f / (1.0f + exp((-60.0f - actualVtilde) / 20.0f));
                real beta_X_r2tilde = 1.12f / (1.0f + exp((actualVtilde - 60.0f) / 20.0f));

                real X_s_inftilde = 1.0f / (1.0f + exp((-5.0f - actualVtilde) / 14.0f));
                real alpha_X_stilde = 1400.0f / sqrt(1.0f + exp((5.0f - actualVtilde) / 6.0f));
                real beta_X_stilde = 1.0f / (1.0f + exp((-35.0f + actualVtilde) / 15.0f));

                real m_inftilde = 1.0f / ((1.0f + exp((-56.86 - actualVtilde) / 9.03f)) * (1.0f + exp((-56.86f - actualVtilde) / 9.03f)));
                real alpha_mtilde = 1.0f / (1.0f + exp((-60.0f - actualVtilde) / 5.0f));
                real beta_mtilde = 0.1f / (1.0f + exp((actualVtilde + 35.0f) / 5.0f)) + (0.1f / (1.0f + exp((actualVtilde - 50.0f) / 200.0f)));

                real h_inftilde = 1.0f / ((1.0f + exp((actualVtilde + 71.55f) / 7.43f)) * (1.0f + exp((actualVtilde + 71.55f) / 7.43f)));
                real alpha_htilde;
                (actualVtilde < -40.0f)
                    ? (alpha_htilde = 0.057f * exp(-(actualVtilde + 80.0f) / 6.8f))
                    : (alpha_htilde = 0.0f);
                real beta_htilde;
                (actualVtilde < -40.0f)
                    ? (beta_htilde = 2.7f * exp(0.079f * actualVtilde) + 3.1f * 1.0e5f * exp(0.3485f * actualVtilde))
                    : (beta_htilde = 0.77f / (0.13f * (1.0f + exp((actualVtilde + 10.66f) / -11.1f))));

                real j_inftilde = 1.0f / ((1.0f + exp((actualVtilde + 71.55f) / 7.43f)) * (1.0f + exp((actualVtilde + 71.55f) / 7.43f)));
                real alpha_jtilde;
                (actualVtilde < -40.0f)
                    ? (alpha_jtilde = ((-25428.0f * exp(0.2444f * actualVtilde) - (6.948e-6f * exp((-0.04391f) * actualVtilde))) * (actualVtilde + 37.78f)) / (1.0f + exp(0.311f * (actualVtilde + 79.23f))))
                    : (alpha_jtilde = 0.0f);
                real beta_jtilde;
                (actualVtilde < -40.0f)
                    ? (beta_jtilde = (0.02424f * exp(-0.01052f * actualVtilde)) / (1.0f + exp(-0.1378f * (actualVtilde + 40.14f))))
                    : (beta_jtilde = (0.6f * exp(0.057f * actualVtilde)) / (1.0f + exp(-0.1f * (actualVtilde + 32.0f))));

                real d_inftilde = 1.0f / (1.0f + exp((-8.0f - actualVtilde) / 7.5f));
                real alpha_dtilde = 1.4f / (1.0f + exp((-35.0f - actualVtilde) / 13.0f)) + 0.25f;
                real beta_dtilde = 1.4f / (1.0f + exp((actualVtilde + 5.0f) / 5.0f));
                real gamma_dtilde = 1.0f / (1.0f + exp((50.0f - actualVtilde) / 20.0f));

                real f_inftilde = 1.0f / (1.0f + exp((actualVtilde + 20.0f) / 7.0f));
                real alpha_ftilde = 1102.5f * exp(-(actualVtilde + 27.0f) * (actualVtilde + 27.0f) / 225.0f);
                real beta_ftilde = 200.0f / (1.0f + exp((13.0f - actualVtilde) / 10.0f));
                real gamma_ftilde = 180.0f / (1.0f + exp((actualVtilde + 30.0f) / 10.0f)) + 20.0f;

                real f2_inftilde = 0.67f / (1.0f + exp((actualVtilde + 35.0f) / 7.0f)) + 0.33f;
                real alpha_f2tilde; // !!!
                alpha_f2tilde = 562.0f * exp(-(actualVtilde + 27.0f) * (actualVtilde + 27.0f) / 240.0f);
                real beta_f2tilde = 31.0f / (1.0f + exp((25.0f - actualVtilde) / 10.0f));
                real gamma_f2tilde; // !!!
                gamma_f2tilde = 80.0f / (1.0f + exp((30.0f + actualVtilde) / 10.0f));

                real fCaSS_inftilde = 0.6f / (1.0f + (Ca_SStilde * Ca_SStilde * 400.0f)) + 0.4f;
                real tau_fCaSStilde = 80.0f / (1.0f + (Ca_SStilde * Ca_SStilde * 400.0f)) + 2.0f;

        #if defined(EPI) || defined(MCELL)
                real s_inftilde = 1.0f / (1.0f + exp((actualVtilde + 20.0f) / 5.0f));
                real tau_stilde = 85.0f * exp(-(actualVtilde + 45.0f) * (actualVtilde + 45.0f) / 320.0f) + 5.0f / (1.0f + exp((actualVtilde - 20.0f) / 5.0f)) + 3.0f;
        #endif // EPI || MCELL
        #ifdef ENDO
                real s_inftilde = 1.0f / (1.0f + exp((actualVtilde + 28.0f) / 5.0f));
                real tau_stilde = 1000.0f * exp(-(actualVtilde + 67.0f) * (actualVtilde + 67.0f) / 1000.0f) + 8.0f;
        #endif // ENDO

                real r_inftilde = 1.0f / (1.0f + exp((20.0f - actualVtilde) / 6.0f));
                real tau_rtilde = 9.5f * exp(-(actualVtilde + 40.0f) * (actualVtilde + 40.0f) / 1800.0f) + 0.8f;

                // Explicit method - auxiliary variables
                real Ileaktilde = V_leak * (Ca_SRtilde - Ca_itilde);
                real Iuptilde = V_maxup / (1.0f + (K_up * K_up) / (Ca_itilde * Ca_itilde));
                real k_CaSRtilde = max_SR - (max_SR - min_SR) / (1.0f + (EC / Ca_SRtilde) * (EC / Ca_SRtilde));
                real k1tilde = k1_prime / k_CaSRtilde;
                real k2tilde = k2_prime * k_CaSRtilde;
                real Otilde = k1tilde * Ca_SStilde * Ca_SStilde * R_primetilde / (k3 + k1tilde * Ca_SStilde * Ca_SStilde);
                real Ireltilde = V_rel * Otilde * (Ca_SRtilde - Ca_SStilde);
                real Ixfertilde = V_xfer * (Ca_SStilde - Ca_itilde);
                real Ca_i_bufCtilde; // !!!
                Ca_i_bufCtilde = 1.0f / (1.0f + bufC * K_bufC / ((K_bufC + Ca_itilde) * (K_bufC + Ca_itilde)));
                real Ca_SR_bufSRtilde; // !!!
                Ca_SR_bufSRtilde = 1.0f / (1.0f + bufSR * K_bufSR / ((K_bufSR + Ca_SRtilde) * (K_bufSR + Ca_SRtilde)));
                real Ca_SS_bufSStilde; // !!!
                Ca_SS_bufSStilde = 1.0f / (1.0f + bufSS * K_bufSS / ((K_bufSS + Ca_SStilde) * (K_bufSS + Ca_SStilde)));

                // Explicit method - RHS of the state variables
                real RHS_R_primetilde_term = (k4 * (1.0f - R_primetilde)) - k2tilde * Ca_SStilde * R_primetilde;
                real RHS_Ca_itilde_term = Ca_i_bufCtilde * ((((Ileaktilde - Iuptilde) * V_SR / V_C) + Ixfertilde) - (((IbCatilde + IpCatilde) - 2.0f * INaCatilde) * Cm / (2.0f * V_C * F)));
                real RHS_Ca_SRtilde_term = Ca_SR_bufSRtilde * (Iuptilde - Ileaktilde - Ireltilde);
                real RHS_Ca_SStilde_term = Ca_SS_bufSStilde * (((-ICaLtilde * Cm / (2.0 * V_SS * F)) + (Ireltilde * V_SR / V_SS)) - (Ixfertilde * V_C / V_SS));
                real RHS_Na_itilde_term = -((INatilde + IbNatilde + 3.0f * INaKtilde + 3.0f * INaCatilde) * Cm / (V_C * F));
                real RHS_K_itilde_term = -((IK1tilde + Itotilde + IKrtilde + IKstilde - 2.0f * INaKtilde + IpKtilde + stim) * Cm / (V_C * F));

                // Rush-Larsen method - update approximations
                X_r1[i] = X_r1_inftilde - (X_r1_inftilde - actualX_r1) * exp(-delta_t / (alpha_X_r1tilde * beta_X_r1tilde));
                X_r2[i] = X_r2_inftilde - (X_r2_inftilde - actualX_r2) * exp(-delta_t / (alpha_X_r2tilde * beta_X_r2tilde));
                X_s[i] = X_s_inftilde - (X_s_inftilde - actualX_s) * exp(-delta_t / (alpha_X_stilde * beta_X_stilde + 80.0f));
                m[i] = m_inftilde - (m_inftilde - actualm) * exp(-delta_t / (alpha_mtilde * beta_mtilde));
                h[i] = h_inftilde - (h_inftilde - actualh) * exp(-delta_t * (alpha_htilde + beta_htilde));
                _j[i] = j_inftilde - (j_inftilde - actualj) * exp(-delta_t * (alpha_jtilde + beta_jtilde));
                d[i] = d_inftilde - (d_inftilde - actuald) * exp(-delta_t / (alpha_dtilde * beta_dtilde + gamma_dtilde));
                f[i] = f_inftilde - (f_inftilde - actualf) * exp(-delta_t / (alpha_ftilde + beta_ftilde + gamma_ftilde));
                f2[i] = f2_inftilde - (f2_inftilde - actualf2) * exp(-delta_t / (alpha_f2tilde + beta_f2tilde + gamma_f2tilde));
                fCaSS[i] = fCaSS_inftilde - (fCaSS_inftilde - actualfCaSS) * exp(-delta_t / tau_fCaSStilde);
                s[i] = s_inftilde - (s_inftilde - actuals) * exp(-delta_t / tau_stilde);
                r[i] = r_inftilde - (r_inftilde - actualr) * exp(-delta_t / tau_rtilde);

                // Explicit method - update approximations
                R_prime[i] = actualR_prime + (delta_t * RHS_R_primetilde_term);
                Ca_i[i] = actualCa_i + (delta_t * RHS_Ca_itilde_term);
                Ca_SR[i] = actualCa_SR + (delta_t * RHS_Ca_SRtilde_term);
                Ca_SS[i] = actualCa_SS + (delta_t * RHS_Ca_SStilde_term);
                Na_i[i] = actualNa_i + (delta_t * RHS_Na_itilde_term);
                K_i[i] = actualK_i + (delta_t * RHS_K_itilde_term);
#endif // TT2
            }

            // ================================================!
            //  Calculate V at n + 1 -> Result goes to V       !
            // ================================================!
            for (int i = 0; i < N; i++)
            {
                real actualV = V[i];
                diff_term = diff_coeff * phi * (V[lim(i - 1, N)] - 2.0f * actualV + V[lim(i + 1, N)]);
                
                LS_b[i] = actualV + (1.0f - theta) * diff_term + partRHS[i];
            }

            tridiag(la, lb, lc, c_prime, d_prime, N, LS_b, result);
            for (int i = 0; i < N; i++)
            {
                V[i] = result[i];
            }

#ifdef SAVE_FRAMES
            // If save frames is true and time step is multiple of frame save rate
            if (timeStepCounter % frameSaveRate == 0)
            {
                // Save frame
                saveFrame(fpFrames, actualTime, V, N);
                printf("Frame at time %.2lf ms saved to %s\n", actualTime, framesPath);
            }
#endif // SAVE_FRAMES

            // Calculate stim velocity
            if (!stim_velocity_measured)
            {
                real begin = L / 3.0f;
                real end = 2.0f * begin;

                if (!aux_stim_velocity_flag)
                {
                    int first_point_index = round(begin / delta_x) + 1;
                    if (V[first_point_index] > 10.0f)
                    {
                        first_point_time = actualTime;
                        aux_stim_velocity_flag = true;
                    }
                }
                else
                {
                    int last_point_index = round(end / delta_x) + 1;
                    if (V[last_point_index] > 10.0f)
                    {
                        last_point_time = actualTime;
                        stim_velocity = (end - begin) / (last_point_time - first_point_time); // cm/ms
                        stim_velocity = stim_velocity * 10.0f; // m/s
                        stim_velocity_measured = true;
                        printf("Stim velocity (measured from %.2f to %.2f cm) is %lf m/s\n", begin, end, stim_velocity);
                    }
                }
            }

            // Update time step counter
            timeStepCounter++;
        }
    }

#endif // not CABLEEQ

    real finishTime = omp_get_wtime();
    real elapsedTime = finishTime - startTime;

// Calculate error
#ifdef CONVERGENCE_ANALYSIS
    real norm2error = calculateNorm2Error(V, exact, N, totalTime, delta_x);
#endif // CONVERGENCE_ANALYSIS

    // Write infos to file
    char infosFilePath[MAX_STRING_SIZE];
    snprintf(infosFilePath, MAX_STRING_SIZE * sizeof(char), "%s/infos/infos_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);
    FILE *fpInfos = fopen(infosFilePath, "w");
    printf("Infos saved to %s\n", infosFilePath);
    fprintf(fpInfos, "Domain Length = %d, Time = %f\n", L, totalTime);
    #ifndef CABLEEQ
    fprintf(fpInfos, "delta_x = %lf, Space steps N = %d, N*N = %d\n", delta_x, N, N * N);
    #else // if not def CABLEEQ
    fprintf(fpInfos, "Stimulus velocity = %lf m/s\n", stim_velocity);
    fprintf(fpInfos, "delta_x = %lf, Space steps N = %d\n", delta_x, N);
    #endif // not CABLEEQ
    fprintf(fpInfos, "delta_t = %lf, Time steps = %d\n", delta_t, M);
    fprintf(fpInfos, "Method %s\n", method);
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
    FILE *fpLastX_r2 = fopen(lastFrameFilePathX_r2, "w");
    printf("Last X_r2 frame saved to %s\n", lastFrameFilePathX_r2);
    FILE *fpLastX_s = fopen(lastFrameFilePathX_s, "w");
    printf("Last X_s frame saved to %s\n", lastFrameFilePathX_s);
    FILE *fpLastm = fopen(lastFrameFilePathm, "w");
    printf("Last m frame saved to %s\n", lastFrameFilePathm);
    FILE *fpLasth = fopen(lastFrameFilePathh, "w");
    printf("Last h frame saved to %s\n", lastFrameFilePathh);
    FILE *fpLastj = fopen(lastFrameFilePathj, "w");
    printf("Last j frame saved to %s\n", lastFrameFilePathj);
    FILE *fpLastd = fopen(lastFrameFilePathd, "w");
    printf("Last d frame saved to %s\n", lastFrameFilePathd);
    FILE *fpLastf = fopen(lastFrameFilePathf, "w");
    printf("Last f frame saved to %s\n", lastFrameFilePathf);
    FILE *fpLastf2 = fopen(lastFrameFilePathf2, "w");
    printf("Last f2 frame saved to %s\n", lastFrameFilePathf2);
    FILE *fpLastfCaSS = fopen(lastFrameFilePathfCaSS, "w");
    printf("Last fCaSS frame saved to %s\n", lastFrameFilePathfCaSS);
    FILE *fpLasts = fopen(lastFrameFilePaths, "w");
    printf("Last s frame saved to %s\n", lastFrameFilePaths);
    FILE *fpLastr = fopen(lastFrameFilePathr, "w");
    printf("Last r frame saved to %s\n", lastFrameFilePathr);
    FILE *fpLastR_prime = fopen(lastFrameFilePathR_prime, "w");
    printf("Last R_prime frame saved to %s\n", lastFrameFilePathR_prime);
    FILE *fpLastCa_i = fopen(lastFrameFilePathCa_i, "w");
    printf("Last Ca_i frame saved to %s\n", lastFrameFilePathCa_i);
    FILE *fpLastCa_SR = fopen(lastFrameFilePathCa_SR, "w");
    printf("Last Ca_SR frame saved to %s\n", lastFrameFilePathCa_SR);
    FILE *fpLastCa_SS = fopen(lastFrameFilePathCa_SS, "w");
    printf("Last Ca_SS frame saved to %s\n", lastFrameFilePathCa_SS);
    FILE *fpLastNa_i = fopen(lastFrameFilePathNa_i, "w");
    printf("Last Na_i frame saved to %s\n", lastFrameFilePathNa_i);
    FILE *fpLastK_i = fopen(lastFrameFilePathK_i, "w");
    printf("Last K_i frame saved to %s\n", lastFrameFilePathK_i);
    #endif // TT2
    #endif // SAVE_LAST_STATE
    
    // Save Action Potential
    #ifdef CABLEEQ
    char APFilePath[MAX_STRING_SIZE];
    snprintf(APFilePath, MAX_STRING_SIZE * sizeof(char), "%s/AP/AP_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);
    FILE *fpAP = fopen(APFilePath, "w");
    printf("Action Potential saved to %s\n", APFilePath);
    #endif // CABLEEQ

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
#ifndef CABLEEQ
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            fprintf(fpLast, "%e ", V[i][j]);
            #ifdef SAVE_LAST_STATE
            #ifdef AFHN
            fprintf(fpLastV, "%e ", V[i][j]);
            fprintf(fpLastW, "%e ", W[i][j]);
            #endif // AFHN
            #ifdef TT2
            fprintf(fpLastV, "%e ", V[i][j]);
            fprintf(fpLastX_r1, "%e ", X_r1[i][j]);
            fprintf(fpLastX_r2, "%e ", X_r2[i][j]);
            fprintf(fpLastX_s, "%e ", X_s[i][j]);
            fprintf(fpLastm, "%e ", m[i][j]);
            fprintf(fpLasth, "%e ", h[i][j]);
            fprintf(fpLastj, "%e ", _j[i][j]);
            fprintf(fpLastd, "%e ", d[i][j]);
            fprintf(fpLastf, "%e ", f[i][j]);
            fprintf(fpLastf2, "%e ", f2[i][j]);
            fprintf(fpLastfCaSS, "%e ", fCaSS[i][j]);
            fprintf(fpLasts, "%e ", s[i][j]);
            fprintf(fpLastr, "%e ", r[i][j]);
            fprintf(fpLastR_prime, "%e ", R_prime[i][j]);
            fprintf(fpLastCa_i, "%e ", Ca_i[i][j]);
            fprintf(fpLastCa_SR, "%e ", Ca_SR[i][j]);
            fprintf(fpLastCa_SS, "%e ", Ca_SS[i][j]);
            fprintf(fpLastNa_i, "%e ", Na_i[i][j]);
            fprintf(fpLastK_i, "%e ", K_i[i][j]);
            #endif // TT2
            #endif // SAVE_LAST_STATE
#ifdef CONVERGENCE_ANALYSIS
            fprintf(fpExact, "%e ", exact[i][j]);
            fprintf(fpErrors, "%e ", abs(V[i][j] - exact[i][j]));
#endif // CONVERGENCE_ANALYSIS
        }
        fprintf(fpLast, "\n");
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
#ifdef CONVERGENCE_ANALYSIS
        fprintf(fpExact, "\n");
        fprintf(fpErrors, "\n");
#endif // CONVERGENCE_ANALYSIS
    }
#else // else of if not CABLEEQ
    for (int i = 0; i < N; i++)
    {
        fprintf(fpLast, "%e ", V[i]);
        #ifdef SAVE_LAST_STATE
        #ifdef AFHN
        fprintf(fpLastV, "%e ", V[i]);
        fprintf(fpLastW, "%e ", W[i]);
        #endif // AFHN
        #ifdef TT2
        fprintf(fpLastV, "%e ", V[i]);
        fprintf(fpLastX_r1, "%e ", X_r1[i]);
        fprintf(fpLastX_r2, "%e ", X_r2[i]);
        fprintf(fpLastX_s, "%e ", X_s[i]);
        fprintf(fpLastm, "%e ", m[i]);
        fprintf(fpLasth, "%e ", h[i]);
        fprintf(fpLastj, "%e ", _j[i]);
        fprintf(fpLastd, "%e ", d[i]);
        fprintf(fpLastf, "%e ", f[i]);
        fprintf(fpLastf2, "%e ", f2[i]);
        fprintf(fpLastfCaSS, "%e ", fCaSS[i]);
        fprintf(fpLasts, "%e ", s[i]);
        fprintf(fpLastr, "%e ", r[i]);
        fprintf(fpLastR_prime, "%e ", R_prime[i]);
        fprintf(fpLastCa_i, "%e ", Ca_i[i]);
        fprintf(fpLastCa_SR, "%e ", Ca_SR[i]);
        fprintf(fpLastCa_SS, "%e ", Ca_SS[i]);
        fprintf(fpLastNa_i, "%e ", Na_i[i]);
        fprintf(fpLastK_i, "%e ", K_i[i]);
        #endif // TT2
        #endif // SAVE_LAST_STATE 
    }

    for (int i = 0; i < M; i++)
    {
        fprintf(fpAP, "%e ", AP[i]);
    }
#endif // not CABLEEQ
    fclose(fpLast);
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
#ifdef CONVERGENCE_ANALYSIS
    fclose(fpExact);
    fclose(fpErrors);
#endif // CONVERGENCE_ANALYSIS

    // Free memory
    free(time);

#ifndef CABLEEQ
    for (int i = 0; i < N; i++)
    {
        free(V[i]);
        free(Vtilde[i]);
        free(RHS[i]);
        free(partRHS[i]);
        free(exact[i]);
#ifdef MONODOMAIN
#ifdef AFHN
        free(W[i]);
#endif // AFHN
#ifdef TT2
        free(X_r1[i]);
        free(X_r2[i]);
        free(X_s[i]);
        free(m[i]);
        free(h[i]);
        free(_j[i]);
        free(d[i]);
        free(f[i]);
        free(f2[i]);
        free(fCaSS[i]);
        free(s[i]);
        free(r[i]);
        free(R_prime[i]);
        free(Ca_i[i]);
        free(Ca_SR[i]);
        free(Ca_SS[i]);
        free(Na_i[i]);
        free(K_i[i]);
#endif // TT2
#endif // MONODOMAIN
    }
#endif // not CABLEEQ
    free(V);
    free(Vtilde);
    free(RHS);
    free(partRHS);
    free(exact);
    free(c_prime);
    free(d_prime);
    free(LS_b);
    free(result);
    free(la);
    free(lb);
    free(lc);
    free(pathToSaveData);
#if defined(MONODOMAIN) || defined(CABLEEQ)
#ifndef CONVERGENCE_ANALYSIS
    free(stimuli);
#endif // not CONVERGENCE_ANALYSIS
#ifdef CABLEEQ
    free(AP);
    fclose(fpAP);
#endif // CABLEEQ
#ifdef AFHN
    free(W);
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
#endif // TT2
#endif // MONODOMAIN || CABLEEQ

    return;
}

#endif // CPU_METHODS_H