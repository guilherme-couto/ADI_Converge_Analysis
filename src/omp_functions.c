#include "../include/core_definitions.h"
#include "../include/config_parser.h"
#include "../include/auxfuncs.h"

#ifdef USE_CUDA

void runSimulationOpenMP(const SimulationConfig *config)
{
    // Unpack configuration parameters
    const ExecutionMode exec_mode = config->exec_mode;
    const CellModel cell_model = config->cell_model;
    const NumericalMethod method = config->method;
    const real delta_t = config->dt;
    const real delta_x = config->dx;
    const real delta_y = config->dy;
    const real totalTime = config->total_time;
    const real Lx = config->Lx;
    const real Ly = config->Ly;
    const int numberOfStimuli = config->stimulus_count;
    const Stimulus *stimuli = config->stimuli;
    const int frameSaveRate = config->frame_save_rate;
    const int numberOfThreads = config->number_of_threads;
    const char *pathToSaveData = config->output_dir;
    const char *pathToRestoreStateFiles = config->path_to_restore_state_files;
    const char *saveFunctionName = config->save_function_name;
    const bool shiftState = config->shift_state;
    const bool saveFrames = config->save_frames;
    const bool saveLastFrame = config->save_last_frame;
    const bool saveLastState = config->save_last_state;
    const bool measureVelocity = config->measure_velocity;

    // Number of steps
    int M = round(totalTime / delta_t); // Number of time steps
    int Nx = round(Lx / delta_x) + 1;   // Spatial steps in x
    printf("Nx = %d\n", Nx);

    int Ny = round(Ly / delta_y) + 1; // Spatial steps in y
    printf("Ny = %d\n", Ny);
    printf("Points in the domain = %d\n", Nx * Ny);

    // Allocate and populate time array
    real *time_array = (real *)malloc(M * sizeof(real));
    initializeTimeArray(time_array, M, delta_t);

    // Allocate 2D arrays for variables
    real **partRHS = (real **)malloc(Ny * sizeof(real *));

#if defined(SSIADI)

    real **RHS = (real **)malloc(Ny * sizeof(real *));

#endif // SSIADI || THETASSIADI

#ifdef AFHN

    real **Vm = (real **)malloc(Ny * sizeof(real *));
    real **W = (real **)malloc(Ny * sizeof(real *));

#endif // AFHN

#ifdef TT2

    real **Vm, **X_r1, **X_r2, **X_s, **m, **h, **_j, **d, **f, **f2, **fCaSS, **s, **r, **Ca_i, **Ca_SR, **Ca_SS, **R_prime, **Na_i, **K_i;
    Vm = (real **)malloc(Ny * sizeof(real *));
    X_r1 = (real **)malloc(Ny * sizeof(real *));
    X_r2 = (real **)malloc(Ny * sizeof(real *));
    X_s = (real **)malloc(Ny * sizeof(real *));
    m = (real **)malloc(Ny * sizeof(real *));
    h = (real **)malloc(Ny * sizeof(real *));
    _j = (real **)malloc(Ny * sizeof(real *));
    d = (real **)malloc(Ny * sizeof(real *));
    f = (real **)malloc(Ny * sizeof(real *));
    f2 = (real **)malloc(Ny * sizeof(real *));
    fCaSS = (real **)malloc(Ny * sizeof(real *));
    s = (real **)malloc(Ny * sizeof(real *));
    r = (real **)malloc(Ny * sizeof(real *));
    Ca_i = (real **)malloc(Ny * sizeof(real *));
    Ca_SR = (real **)malloc(Ny * sizeof(real *));
    Ca_SS = (real **)malloc(Ny * sizeof(real *));
    R_prime = (real **)malloc(Ny * sizeof(real *));
    Na_i = (real **)malloc(Ny * sizeof(real *));
    K_i = (real **)malloc(Ny * sizeof(real *));

#endif // TT2

#ifdef MV

    real **Vm, **v, **w, **s;
    Vm = (real **)malloc(Ny * sizeof(real *));
    v = (real **)malloc(Ny * sizeof(real *));
    w = (real **)malloc(Ny * sizeof(real *));
    s = (real **)malloc(Ny * sizeof(real *));

#endif // MV

    for (int i = 0; i < Ny; i++)
    {
        partRHS[i] = (real *)malloc(Nx * sizeof(real));

#if defined(SSIADI)

        RHS[i] = (real *)malloc(Nx * sizeof(real));

#endif // SSIADI || THETASSIADI

#ifdef AFHN

        Vm[i] = (real *)malloc(Nx * sizeof(real));
        W[i] = (real *)malloc(Nx * sizeof(real));

#endif // AFHN

#ifdef TT2

        Vm[i] = (real *)malloc(Nx * sizeof(real));
        X_r1[i] = (real *)malloc(Nx * sizeof(real));
        X_r2[i] = (real *)malloc(Nx * sizeof(real));
        X_s[i] = (real *)malloc(Nx * sizeof(real));
        m[i] = (real *)malloc(Nx * sizeof(real));
        h[i] = (real *)malloc(Nx * sizeof(real));
        _j[i] = (real *)malloc(Nx * sizeof(real));
        d[i] = (real *)malloc(Nx * sizeof(real));
        f[i] = (real *)malloc(Nx * sizeof(real));
        f2[i] = (real *)malloc(Nx * sizeof(real));
        fCaSS[i] = (real *)malloc(Nx * sizeof(real));
        s[i] = (real *)malloc(Nx * sizeof(real));
        r[i] = (real *)malloc(Nx * sizeof(real));
        Ca_i[i] = (real *)malloc(Nx * sizeof(real));
        Ca_SR[i] = (real *)malloc(Nx * sizeof(real));
        Ca_SS[i] = (real *)malloc(Nx * sizeof(real));
        R_prime[i] = (real *)malloc(Nx * sizeof(real));
        Na_i[i] = (real *)malloc(Nx * sizeof(real));
        K_i[i] = (real *)malloc(Nx * sizeof(real));

#endif // TT2

#ifdef MV

        Vm[i] = (real *)malloc(Nx * sizeof(real));
        v[i] = (real *)malloc(Nx * sizeof(real));
        w[i] = (real *)malloc(Nx * sizeof(real));
        s[i] = (real *)malloc(Nx * sizeof(real));

#endif // MV
    }

    // Initialize variables
#ifdef AFHN

    initialize2DVariableWithValue(Vm, Nx, Ny, Vm_init);
    initialize2DVariableWithValue(W, Nx, Ny, W_init);

#endif // AFHN

#ifdef TT2

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

    initialize2DVariableWithValue(Vm, Nx, Ny, u_init);
    initialize2DVariableWithValue(v, Nx, Ny, v_init);
    initialize2DVariableWithValue(w, Nx, Ny, w_init);
    initialize2DVariableWithValue(s, Nx, Ny, s_init);

#endif // MV

#ifdef RESTORE_STATE

    printf("\n");
    printf("Restoring state variables...\n");

    // Initialize variables with a solution
    real real_ref_dx = 0.002f;
    real real_def_dy = 0.002f;

    char *pathToRestoreStateFiles = (char *)malloc(MAX_STRING_SIZE);

#ifdef AFHN

    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeVm.txt", REAL_TYPE, PROBLEM, cell_model);
    initialize2DVariableFromFile(Vm, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "Vm", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeW.txt", REAL_TYPE, PROBLEM, cell_model);
    initialize2DVariableFromFile(W, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "W", real_ref_dx, real_def_dy);

#endif // AFHN

#ifdef TT2

    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeVm.txt", REAL_TYPE, PROBLEM, cell_model);
    initialize2DVariableFromFile(Vm, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "Vm", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeX_r1.txt", REAL_TYPE, PROBLEM, cell_model);
    initialize2DVariableFromFile(X_r1, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "X_r1", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeX_r2.txt", REAL_TYPE, PROBLEM, cell_model);
    initialize2DVariableFromFile(X_r2, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "X_r2", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeX_s.txt", REAL_TYPE, PROBLEM, cell_model);
    initialize2DVariableFromFile(X_s, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "X_s", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframem.txt", REAL_TYPE, PROBLEM, cell_model);
    initialize2DVariableFromFile(m, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "m", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeh.txt", REAL_TYPE, PROBLEM, cell_model);
    initialize2DVariableFromFile(h, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "h", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframej.txt", REAL_TYPE, PROBLEM, cell_model);
    initialize2DVariableFromFile(_j, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "j", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframed.txt", REAL_TYPE, PROBLEM, cell_model);
    initialize2DVariableFromFile(d, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "d", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframef.txt", REAL_TYPE, PROBLEM, cell_model);
    initialize2DVariableFromFile(f, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "f", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframef2.txt", REAL_TYPE, PROBLEM, cell_model);
    initialize2DVariableFromFile(f2, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "f2", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframefCaSS.txt", REAL_TYPE, PROBLEM, cell_model);
    initialize2DVariableFromFile(fCaSS, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "fCaSS", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframes.txt", REAL_TYPE, PROBLEM, cell_model);
    initialize2DVariableFromFile(s, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "s", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframer.txt", REAL_TYPE, PROBLEM, cell_model);
    initialize2DVariableFromFile(r, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "r", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeCa_i.txt", REAL_TYPE, PROBLEM, cell_model);
    initialize2DVariableFromFile(Ca_i, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "Ca_i", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeCa_SR.txt", REAL_TYPE, PROBLEM, cell_model);
    initialize2DVariableFromFile(Ca_SR, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "Ca_SR", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeCa_SS.txt", REAL_TYPE, PROBLEM, cell_model);
    initialize2DVariableFromFile(Ca_SS, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "Ca_SS", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeR_prime.txt", REAL_TYPE, PROBLEM, cell_model);
    initialize2DVariableFromFile(R_prime, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "R_prime", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeNa_i.txt", REAL_TYPE, PROBLEM, cell_model);
    initialize2DVariableFromFile(Na_i, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "Na_i", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeK_i.txt", REAL_TYPE, PROBLEM, cell_model);
    initialize2DVariableFromFile(K_i, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "K_i", real_ref_dx, real_def_dy);

#endif // TT2

#ifdef MV

    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframeVm.txt", REAL_TYPE, PROBLEM, cell_model);
    initialize2DVariableFromFile(Vm, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "Vm", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframev.txt", REAL_TYPE, PROBLEM, cell_model);
    initialize2DVariableFromFile(v, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "v", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframew.txt", REAL_TYPE, PROBLEM, cell_model);
    initialize2DVariableFromFile(w, Nx, Ny, pathToRestoreStateFiles, delta_x, delta_y, "w", real_ref_dx, real_def_dy);
    snprintf(pathToRestoreStateFiles, MAX_STRING_SIZE, "./restore_state/%s/%s/%s/lastframes.txt", REAL_TYPE, PROBLEM, cell_model);
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
    real diff_coeff;

#ifdef AFHN

    diff_coeff = sigma / (Cm * chi);

#endif // AFHN

#ifdef TT2

    diff_coeff = sigma / chi;

#endif // TT2

#ifdef MV

    diff_coeff = Dtilde;

#endif // MV

#if defined(SSIADI) || defined(OSADI)

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

#ifdef OSADI

    populateDiagonalThomasAlgorithm(la_x, lb_x, lc_x, Nx, phi_x * diff_coeff);
    populateDiagonalThomasAlgorithm(la_y, lb_y, lc_y, Ny, phi_y * diff_coeff);

#endif // OSADI

#endif // SSIADI || OSADI

#if defined(MONODOMAIN)

#endif // MONODOMAIN

    // Create directories
    char *pathToSaveData = (char *)malloc(MAX_STRING_SIZE);
    createDirectories(delta_t, delta_x, delta_y, pathToSaveData);

#ifdef SAVE_FRAMES

    // Save frames variables
    char framesPath[MAX_STRING_SIZE];
    FILE *fpFrames;
    snprintf(framesPath, MAX_STRING_SIZE, "%s/frames.txt", pathToSaveData);
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
    real startTime, finishTime, elapsedTime1stPart, elapsedTime2ndPart = 0.0f;
    real startLSTime, finishLSTime, elapsedTime1stLS, elapsedTime2ndLS = 0.0f;

    // Auxiliary variables for the time loop
    int timeStepCounter = 0;
    real actualTime = 0.0f;

    omp_set_num_threads(NUMTHREADS);
    printf("Number of threads = %d\n", NUMTHREADS);

    // Private variables for OpenMP
    int i, j, si;
    real actualVm;
    real diff_term = 0.0f;
    real stim = 0.0f;
    real point_potential;

    // Start measuring the execution time
    printf("\n");
    printf("Starting simulation...\n");
    startExecutionTime = omp_get_wtime();

#ifdef AFHN

    real actualW, RHS_Vm_term;
    real Vmtilde, Wtilde, RHS_Vmtilde_term;

#pragma omp parallel default(shared) \
    private(i, j, si, actualVm, diff_term, stim, actualW, RHS_Vm_term, Vmtilde, Wtilde, RHS_Vmtilde_term)

#endif // AFHN

#ifdef TT2

    // TODO: Implement

#endif // TT2

#ifdef MV

    real actualv, actualw, actuals;
    real Htheta_w, Htheta_o, Htheta_v, Htheta_vminus;
    real tau_o, tau_so, tau_vminus;
    real J_fi, J_so, J_si;
    real RHS_Vm_term;
    real Vmtilde, vtilde, wtilde, stilde;
    real v_inf, tau_v_RL, v_inf_RL;
    real w_inf, tau_wminus, tau_w_RL, w_inf_RL;
    real tau_s, s_inf_RL;
    real J_fi_tilde, J_so_tilde, J_si_tilde;
    real RHS_Vmtilde_term;

#pragma omp parallel default(shared)                                                                              \
    private(i, j, si, actualVm, diff_term, stim, actualv, actualw, actuals, Htheta_w, Htheta_o, Htheta_v,         \
                Htheta_vminus, tau_o, tau_so, tau_vminus, J_fi, J_so, J_si, RHS_Vm_term, Vmtilde, vtilde, wtilde, \
                stilde, v_inf, tau_v_RL, v_inf_RL, w_inf, tau_wminus, tau_w_RL, w_inf_RL, tau_s, s_inf_RL)

#endif // MV

    while (timeStepCounter < M)
    {
        // Get time step
        actualTime = time_array[timeStepCounter];

// Start measuring time of 1st part
#pragma omp barrier
#pragma omp single
        startTime = omp_get_wtime();

// ================================================!
//  Calculate Approxs. and Update ODEs             !
// ================================================!
#pragma omp for collapse(2) schedule(static)
        for (i = 0; i < Ny; i++)
        {
            for (j = 0; j < Nx; j++)
            {
#ifdef AFHN

                // Calculate the explicit part of the RHS, including the diffusion term in both directions
                actualVm = Vm[i][j];
                actualW = W[i][j];
                RHS_Vm_term = (G * actualVm * (1.0f - (actualVm / vth)) * (1.0f - (actualVm / vp))) + (eta1 * actualVm * actualW);

                // Stimulation
                stim = 0.0f;

#pragma unroll
                for (si = 0; si < numberOfStimuli; si++)
                {
                    if (actualTime >= stimuli[si].begin_time && actualTime <= stimuli[si].begin_time + stimuli[si].duration && j >= stimuli[si].xMinDisc && j <= stimuli[si].xMaxDisc && i >= stimuli[si].yMinDisc && i <= stimuli[si].yMaxDisc)
                    {
                        stim = stimuli[si].amplitude;
                        break;
                    }
                }

#if defined(SSIADI)

                // Calculate aproximation with RK2 -> Vmn+1/2 = Vmn + 0.5*diffusion + 0.5*dt*R(Vmn, Wn)
                diff_term = diff_coeff * (phi_x * (Vm[i][lim(j - 1, Nx)] - 2.0f * actualVm + Vm[i][lim(j + 1, Nx)]) + phi_y * (Vm[lim(i - 1, Ny)][j] - 2.0f * actualVm + Vm[lim(i + 1, Ny)][j]));
                Vmtilde = actualVm + 0.5f * diff_term + (0.5f * delta_t * (stim - RHS_Vm_term));

                // Calculate approximation for state variables and prepare part of the RHS of the following linear systems
                Wtilde = actualW + (0.5f * delta_t * (eta2 * ((actualVm / vp) - (eta3 * actualW))));
                RHS_Vmtilde_term = (G * Vmtilde * (1.0f - (Vmtilde / vth)) * (1.0f - (Vmtilde / vp))) + (eta1 * Vmtilde * Wtilde);
                partRHS[i][j] = delta_t * (stim - RHS_Vmtilde_term);

                // Update state variables
                W[i][j] = actualW + delta_t * (eta2 * ((Vmtilde / vp) - (eta3 * Wtilde))); // with RK2 -> Wn+1 = Wn + dt*R(Vm*, W*)

#endif // SSIADI

#if defined(OSADI)

                // Calculate part of the RHS of the following linear systems with Forward Euler
                partRHS[i][j] = delta_t * (stim - RHS_Vm_term);

                // Update state variables
                W[i][j] = actualW + delta_t * (eta2 * ((actualVm / vp) - (eta3 * actualW))); // with Forward Euler -> Wn+1 = Wn + dt*R(Vmn, Wn)

#endif // OSADI

#ifdef FE

                // Update variables explicitly
                diff_term = diff_coeff * (phi_x * (Vm[i][lim(j - 1, Nx)] - 2.0f * actualVm + Vm[i][lim(j + 1, Nx)]) + phi_y * (Vm[lim(i - 1, Ny)][j] - 2.0f * actualVm + Vm[lim(i + 1, Ny)][j]));
                partRHS[i][j] = actualVm + diff_term + delta_t * (stim - RHS_Vm_term);

                W[i][j] = actualW + delta_t * (eta2 * ((actualVm / vp) - (eta3 * actualW)));

#endif // FE

#endif // AFHN

#ifdef TT2
                // TODO: Implement TT2 for MONODOMAIN in SERIAL mode
#endif // TT2

#ifdef MV

                // Calculate the explicit part of the RHS, including the diffusion term in both directions
                actualVm = Vm[i][j];
                actualv = v[i][j];
                actualw = w[i][j];
                actuals = s[i][j];

                stim = 0.0f;

#pragma unroll
                for (si = 0; si < numberOfStimuli; si++)
                {
                    if (actualTime >= stimuli[si].begin_time && actualTime <= stimuli[si].begin_time + stimuli[si].duration && j >= stimuli[si].xMinDisc && j <= stimuli[si].xMaxDisc && i >= stimuli[si].yMinDisc && i <= stimuli[si].yMaxDisc)
                    {
                        stim = stimuli[si].amplitude;
                        break;
                    }
                }

                diff_term = diff_coeff * (phi_x * (Vm[i][lim(j - 1, Nx)] - 2.0f * actualVm + Vm[i][lim(j + 1, Nx)]) + phi_y * (Vm[lim(i - 1, Ny)][j] - 2.0f * actualVm + Vm[lim(i + 1, Ny)][j]));

                // Calculate RHS of the equations
                // Auxiliary variables
                Htheta_w = (actualVm - theta_w > 0.0f) ? 1.0f : 0.0f;
                Htheta_o = (actualVm - theta_o > 0.0f) ? 1.0f : 0.0f;
                Htheta_v = (actualVm - theta_v > 0.0f) ? 1.0f : 0.0f;
                Htheta_vminus = (actualVm - theta_vminus > 0.0f) ? 1.0f : 0.0f;

                tau_o = (1.0f - Htheta_o) * tau_o1 + Htheta_o * tau_o2;
                tau_so = tau_so1 + (tau_so2 - tau_so1) * (1.0f + tanh(k_so * (actualVm - u_so))) * 0.5f;
                tau_vminus = (1.0f - Htheta_vminus) * tau_v1minus + Htheta_vminus * tau_v2minus;

                // Currents
                J_fi = -actualv * Htheta_v * (actualVm - theta_v) * (u_u - actualVm) / tau_fi;
                J_so = ((actualVm - u_o) * (1.0f - Htheta_w) / tau_o) + (Htheta_w / tau_so);
                J_si = -Htheta_w * actualw * actuals / tau_si;

                // RHS of the state variables
                RHS_Vm_term = J_fi + J_so + J_si;

#if defined(SSIADI)

                // Calculate Vmtilde -> utilde = u^n + 0.5 * dt * (A*u^n + R(u^n))
                Vmtilde = actualVm + 0.5f * diff_term + (0.5f * delta_t * (stim - RHS_Vm_term));

                // Auxiliary variables for Rush-Larsen or Forward Euler
                v_inf = ((actualVm < theta_vminus) ? 1.0f : 0.0f);
                tau_v_RL = (tau_vplus * tau_vminus) / (tau_vplus - tau_vplus * Htheta_v + tau_vminus * Htheta_v);
                v_inf_RL = (tau_vplus * v_inf * (1.0f - Htheta_v)) / (tau_vplus - tau_vplus * Htheta_v + tau_vminus * Htheta_v);

                w_inf = ((1.0f - Htheta_o) * (1.0f - (actualVm / tau_winf)) + Htheta_o * w_infstar);
                tau_wminus = tau_w1minus + (tau_w2minus - tau_w1minus) * (1.0f + tanh(k_wminus * (actualVm - u_wminus))) * 0.5f;
                tau_w_RL = (tau_wplus * tau_wminus) / (tau_wplus - tau_wplus * Htheta_w + tau_wminus * Htheta_w);
                w_inf_RL = (tau_wplus * w_inf * (1.0f - Htheta_w)) / (tau_wplus - tau_wplus * Htheta_w + tau_wminus * Htheta_w);

                tau_s = (1.0f - Htheta_w) * tau_s1 + Htheta_w * tau_s2;
                s_inf_RL = (1.0f + tanh(k_s * (actualVm - u_s))) * 0.5f;

                // Calculate approximations with Rush-Larsen or Forward Euler using half time step
                vtilde, wtilde, stilde;
                (tau_v_RL > 1.0e-10)
                    ? (vtilde = v_inf_RL - (v_inf_RL - actualv) * exp(-0.5f * delta_t / tau_v_RL))
                    : (vtilde = actualv + 0.5f * delta_t * (1.0f - Htheta_v) * (v_inf - actualv) / tau_vminus - Htheta_v * actualv / tau_vplus);

                (tau_w_RL > 1.0e-10)
                    ? (wtilde = w_inf_RL - (w_inf_RL - actualw) * exp(-0.5f * delta_t / tau_w_RL))
                    : (wtilde = actualw + 0.5f * delta_t * (1.0f - Htheta_w) * (w_inf - actualw) / tau_wminus - Htheta_w * actualw / tau_wplus);

                (tau_s > 1.0e-10)
                    ? (stilde = s_inf_RL - (s_inf_RL - actuals) * exp(-0.5f * delta_t / tau_s))
                    : (stilde = actuals + 0.5f * delta_t * (s_inf_RL - actuals) / tau_s);

                // Calculate RHS of the equations with approximations
                // Auxiliary variables
                Htheta_w = (Vmtilde - theta_w > 0.0f) ? 1.0f : 0.0f;
                Htheta_o = (Vmtilde - theta_o > 0.0f) ? 1.0f : 0.0f;
                Htheta_v = (Vmtilde - theta_v > 0.0f) ? 1.0f : 0.0f;
                Htheta_vminus = (Vmtilde - theta_vminus > 0.0f) ? 1.0f : 0.0f;

                tau_o = (1.0f - Htheta_o) * tau_o1 + Htheta_o * tau_o2;
                tau_so = tau_so1 + (tau_so2 - tau_so1) * (1.0f + tanh(k_so * (Vmtilde - u_so))) * 0.5f;
                tau_vminus = (1.0f - Htheta_vminus) * tau_v1minus + Htheta_vminus * tau_v2minus;

                // Currents
                J_fi_tilde = -actualv * Htheta_v * (Vmtilde - theta_v) * (u_u - Vmtilde) / tau_fi;
                J_so_tilde = ((Vmtilde - u_o) * (1.0f - Htheta_w) / ((1.0f - Htheta_o) * tau_o1 + Htheta_o * tau_o2)) + (Htheta_w / (tau_so1 + (((tau_so2 - tau_so1) * (1.0f + tanh(k_so * (Vmtilde - u_so)))) * 0.5f)));
                J_si_tilde = -Htheta_w * wtilde * stilde / tau_si;

                // Update partRHS
                RHS_Vmtilde_term = J_fi_tilde + J_so_tilde + J_si_tilde;
                partRHS[i][j] = delta_t * (stim - RHS_Vmtilde_term);

                // Update auxiliary variables for Rush-Larsen or Forward Euler with approximations
                v_inf = ((Vmtilde < theta_vminus) ? 1.0f : 0.0f);
                tau_v_RL = (tau_vplus * tau_vminus) / (tau_vplus - tau_vplus * Htheta_v + tau_vminus * Htheta_v);
                v_inf_RL = (tau_vplus * v_inf * (1.0f - Htheta_v)) / (tau_vplus - tau_vplus * Htheta_v + tau_vminus * Htheta_v);

                w_inf = ((1.0f - Htheta_o) * (1.0f - (Vmtilde / tau_winf)) + Htheta_o * w_infstar);
                tau_wminus = tau_w1minus + (tau_w2minus - tau_w1minus) * (1.0f + tanh(k_wminus * (Vmtilde - u_wminus))) * 0.5f;
                tau_w_RL = (tau_wplus * tau_wminus) / (tau_wplus - tau_wplus * Htheta_w + tau_wminus * Htheta_w);
                w_inf_RL = (tau_wplus * w_inf * (1.0f - Htheta_w)) / (tau_wplus - tau_wplus * Htheta_w + tau_wminus * Htheta_w);

                tau_s = (1.0f - Htheta_w) * tau_s1 + Htheta_w * tau_s2;
                s_inf_RL = (1.0f + tanh(k_s * (Vmtilde - u_s))) * 0.5f;

                // Update state variables with Rush-Larsen or RK2 (Heun's method) using approximations
                (tau_v_RL > 1.0e-10)
                    ? (v[i][j] = v_inf_RL - (v_inf_RL - actualv) * exp(-delta_t / tau_v_RL))
                    : (v[i][j] = actualv + delta_t * (1.0f - Htheta_v) * (v_inf - vtilde) / tau_vminus - Htheta_v * vtilde / tau_vplus);

                (tau_w_RL > 1.0e-10)
                    ? (w[i][j] = w_inf_RL - (w_inf_RL - actualw) * exp(-delta_t / tau_w_RL))
                    : (w[i][j] = actualw + delta_t * (1.0f - Htheta_w) * (w_inf - wtilde) / tau_wminus - Htheta_w * wtilde / tau_wplus);

                (tau_s > 1.0e-10)
                    ? (s[i][j] = s_inf_RL - (s_inf_RL - actuals) * exp(-delta_t / tau_s))
                    : (s[i][j] = actuals + delta_t * (s_inf_RL - stilde) / tau_s);

#endif // SSIADI

#ifdef OSADI

                // Update partRHS
                partRHS[i][j] = delta_t * (stim - RHS_Vm_term);

                // Auxiliary variables for Rush-Larsen or Forward Euler
                v_inf = ((actualVm < theta_vminus) ? 1.0f : 0.0f);
                tau_v_RL = (tau_vplus * tau_vminus) / (tau_vplus - tau_vplus * Htheta_v + tau_vminus * Htheta_v);
                v_inf_RL = (tau_vplus * v_inf * (1.0f - Htheta_v)) / (tau_vplus - tau_vplus * Htheta_v + tau_vminus * Htheta_v);

                w_inf = ((1.0f - Htheta_o) * (1.0f - (actualVm / tau_winf)) + Htheta_o * w_infstar);
                tau_wminus = tau_w1minus + (tau_w2minus - tau_w1minus) * (1.0f + tanh(k_wminus * (actualVm - u_wminus))) * 0.5f;
                tau_w_RL = (tau_wplus * tau_wminus) / (tau_wplus - tau_wplus * Htheta_w + tau_wminus * Htheta_w);
                w_inf_RL = (tau_wplus * w_inf * (1.0f - Htheta_w)) / (tau_wplus - tau_wplus * Htheta_w + tau_wminus * Htheta_w);

                tau_s = (1.0f - Htheta_w) * tau_s1 + Htheta_w * tau_s2;
                s_inf_RL = (1.0f + tanh(k_s * (actualVm - u_s))) * 0.5f;

                // Update state variables with Rush-Larsen or Forward Euler
                (tau_v_RL > 1.0e-10)
                    ? (v[i][j] = v_inf_RL - (v_inf_RL - actualv) * exp(-delta_t / tau_v_RL))
                    : (v[i][j] = actualv + delta_t * (1.0f - Htheta_v) * (v_inf - actualv) / tau_vminus - Htheta_v * actualv / tau_vplus);

                (tau_w_RL > 1.0e-10)
                    ? (w[i][j] = w_inf_RL - (w_inf_RL - actualw) * exp(-delta_t / tau_w_RL))
                    : (w[i][j] = actualw + delta_t * (1.0f - Htheta_w) * (w_inf - actualw) / tau_wminus - Htheta_w * actualw / tau_wplus);

                (tau_s > 1.0e-10)
                    ? (s[i][j] = s_inf_RL - (s_inf_RL - actuals) * exp(-delta_t / tau_s))
                    : (s[i][j] = actuals + delta_t * (s_inf_RL - actuals) / tau_s);

#endif // OSADI

#ifdef FE

                // Update partRHS (auxiliary variable) with Forward Euler
                diff_term = diff_coeff * (phi_x * (Vm[i][lim(j - 1, Nx)] - 2.0f * actualVm + Vm[i][lim(j + 1, Nx)]) + phi_y * (Vm[lim(i - 1, Ny)][j] - 2.0f * actualVm + Vm[lim(i + 1, Ny)][j]));
                partRHS[i][j] = actualVm + diff_term + delta_t * (stim - RHS_Vm_term);

                // Auxiliary variables for Rush-Larsen or Forward Euler
                v_inf = ((actualVm < theta_vminus) ? 1.0f : 0.0f);
                tau_v_RL = (tau_vplus * tau_vminus) / (tau_vplus - tau_vplus * Htheta_v + tau_vminus * Htheta_v);
                v_inf_RL = (tau_vplus * v_inf * (1.0f - Htheta_v)) / (tau_vplus - tau_vplus * Htheta_v + tau_vminus * Htheta_v);

                w_inf = ((1.0f - Htheta_o) * (1.0f - (actualVm / tau_winf)) + Htheta_o * w_infstar);
                tau_wminus = tau_w1minus + (tau_w2minus - tau_w1minus) * (1.0f + tanh(k_wminus * (actualVm - u_wminus))) * 0.5f;
                tau_w_RL = (tau_wplus * tau_wminus) / (tau_wplus - tau_wplus * Htheta_w + tau_wminus * Htheta_w);
                w_inf_RL = (tau_wplus * w_inf * (1.0f - Htheta_w)) / (tau_wplus - tau_wplus * Htheta_w + tau_wminus * Htheta_w);

                tau_s = (1.0f - Htheta_w) * tau_s1 + Htheta_w * tau_s2;
                s_inf_RL = (1.0f + tanh(k_s * (actualVm - u_s))) * 0.5f;

                // Update state variables with Rush-Larsen or Forward Euler
                (tau_v_RL > 1.0e-10)
                    ? (v[i][j] = v_inf_RL + (actualv - v_inf_RL) * exp(-delta_t / tau_v_RL))
                    : (v[i][j] = actualv + delta_t * ((1.0f - Htheta_v) * (v_inf - actualv) / tau_vminus - Htheta_v * actualv / tau_vplus));

                (tau_w_RL > 1.0e-10)
                    ? (w[i][j] = w_inf_RL + (actualw - w_inf_RL) * exp(-delta_t / tau_w_RL))
                    : (w[i][j] = actualw + delta_t * ((1.0f - Htheta_w) * (w_inf - actualw) / tau_wminus - Htheta_w * actualw / tau_wplus));

                (tau_s > 1.0e-10)
                    ? (s[i][j] = s_inf_RL + (actuals - s_inf_RL) * exp(-delta_t / tau_s))
                    : (s[i][j] = actuals + delta_t * (s_inf_RL - actuals) / tau_s);

#endif // FE

#endif // MV
            }
        }

#pragma omp barrier
#pragma omp single
        {
            // End of the 1st part of the time step
            finishTime = omp_get_wtime();
            elapsedTime1stPart += finishTime - startTime;

            // Start measuring time of 2nd part
            startTime = omp_get_wtime();
        }

#if defined(SSIADI)

        real *c_prime_y = (real *)malloc(Ny * sizeof(real));
        real *d_prime_y = (real *)malloc(Ny * sizeof(real));
        real *LS_b_y = (real *)malloc(Ny * sizeof(real));
        real *result_y = (real *)malloc(Ny * sizeof(real));

// ================================================!
//  Calculate Vm at n+1/2 -> Result goes to RHS    !
//  diffusion implicit in y and explicit in x      !
// ================================================!
#pragma omp barrier
#pragma omp for schedule(static)
        for (j = 0; j < Nx; j++)
        {
            // Calculate the RHS of the linear system with the explicit diffusion term along x
            for (i = 0; i < Ny; i++)
            {
                actualVm = Vm[i][j];
                diff_term = diff_coeff * tau * phi_x * (Vm[i][lim(j - 1, Nx)] - 2.0f * actualVm + Vm[i][lim(j + 1, Nx)]);
                LS_b_y[i] = actualVm + diff_term + 0.5f * partRHS[i][j];
            }

            // Solve the linear system
            tridiagonalSystemSolver(la_y, lb_y, lc_y, c_prime_y, d_prime_y, Ny, LS_b_y, result_y);

            // Update with the result
            for (i = 0; i < Ny; i++)
                RHS[i][j] = result_y[i];
        }

        free(c_prime_y);
        free(d_prime_y);
        free(LS_b_y);
        free(result_y);

        real *c_prime_x = (real *)malloc(Nx * sizeof(real));
        real *d_prime_x = (real *)malloc(Nx * sizeof(real));
        real *LS_b_x = (real *)malloc(Nx * sizeof(real));
        real *result_x = (real *)malloc(Nx * sizeof(real));

// ================================================!
//  Calculate Vm at n+1 -> Result goes to Vm       !
//  diffusion implicit in x and explicit in y      !
// ================================================!
#pragma omp barrier
#pragma omp for schedule(static)
        for (i = 0; i < Ny; i++)
        {
            // Calculate the RHS of the linear system with the explicit diffusion term along y
            for (j = 0; j < Nx; j++)
            {
                actualVm = RHS[i][j];
                diff_term = diff_coeff * tau * phi_y * (RHS[lim(i - 1, Ny)][j] - 2.0f * actualVm + RHS[lim(i + 1, Ny)][j]);
                LS_b_x[j] = actualVm + diff_term + 0.5f * partRHS[i][j];
            }

            // Solve the linear system
            tridiagonalSystemSolver(la_x, lb_x, lc_x, c_prime_x, d_prime_x, Nx, LS_b_x, result_x);

            // Update with the result
            for (j = 0; j < Nx; j++)
                Vm[i][j] = result_x[j];
        }

        free(c_prime_x);
        free(d_prime_x);
        free(LS_b_x);
        free(result_x);

#endif // SSIADI

#ifdef OSADI

        real *c_prime_y = (real *)malloc(Ny * sizeof(real));
        real *d_prime_y = (real *)malloc(Ny * sizeof(real));
        real *LS_b_y = (real *)malloc(Ny * sizeof(real));
        real *result_y = (real *)malloc(Ny * sizeof(real));

// ================================================!
//  Calculate Vm at n+1/2 -> Result goes to Vm      !
// ================================================!
#pragma omp barrier
#pragma omp for schedule(static)
        for (j = 0; j < Nx; j++)
        {
            // Calculate the RHS of the linear system
            for (i = 0; i < Ny; i++)
                LS_b_y[i] = Vm[i][j] + 0.5f * partRHS[i][j];

            // Solve the linear system
            tridiagonalSystemSolver(la_y, lb_y, lc_y, c_prime_y, d_prime_y, Ny, LS_b_y, result_y);

            // Update with the result
            for (i = 0; i < Ny; i++)
                Vm[i][j] = result_y[i];
        }

        free(c_prime_y);
        free(d_prime_y);
        free(LS_b_y);
        free(result_y);

        real *c_prime_x = (real *)malloc(Nx * sizeof(real));
        real *d_prime_x = (real *)malloc(Nx * sizeof(real));
        real *LS_b_x = (real *)malloc(Nx * sizeof(real));
        real *result_x = (real *)malloc(Nx * sizeof(real));

// ================================================!
//  Calculate Vm at n+1 -> Result goes to Vm       !
// ================================================!
#pragma omp barrier
#pragma omp for schedule(static)
        for (i = 0; i < Ny; i++)
        {
            // Calculate the RHS of the linear system
            for (j = 0; j < Nx; j++)
                LS_b_x[j] = Vm[i][j] + 0.5f * partRHS[i][j];

            // Solve the linear system
            tridiagonalSystemSolver(la_x, lb_x, lc_x, c_prime_x, d_prime_x, Nx, LS_b_x, result_x);

            // Update with the result
            for (j = 0; j < Nx; j++)
                Vm[i][j] = result_x[j];
        }

        free(c_prime_x);
        free(d_prime_x);
        free(LS_b_x);
        free(result_x);

#endif // OSADI

#ifdef FE

// ==================!
//  Update Vm        !
// ==================!
#pragma omp barrier
#pragma omp for collapse(2) schedule(static)
        for (i = 0; i < Ny; i++)
            for (j = 0; j < Nx; j++)
                Vm[i][j] = partRHS[i][j];

#endif // FE

// End of the 2nd part of the time step
#pragma omp barrier
#pragma omp single
        {
            finishTime = omp_get_wtime();
            elapsedTime2ndPart += finishTime - startTime;
        }

#ifdef SAVE_FRAMES

#pragma omp single
        {
            startSaveFramesTime = omp_get_wtime();

            // If save frames is true and time step is multiple of frame save rate
            if (timeStepCounter % frameSaveRate == 0)
            {
                // Save frame
                fprintf(fpFrames, "%lf\n", actualTime);
                saveFrame(fpFrames, Vm, Nx, Ny);
                SUCCESSMSG("Frame at time %.2f ms saved to %s\n", actualTime, framesPath);
            }

            finishSaveFramesTime = omp_get_wtime();
            elapsedSaveFramesTime += finishSaveFramesTime - startSaveFramesTime;
        }

#endif // SAVE_FRAMES

#ifdef MEASURE_VELOCITY

#pragma omp single
        {
            startMeasureVelocityTime = omp_get_wtime();

            // Calculate stim velocity
            if (!stim_velocity_measured)
            {
                real point_potential = 0.0f;
                if (!aux_stim_velocity_flag)
                {
                    point_potential = Vm[0][begin_point_index];

#ifdef MV

                    point_potential = rescaleVm(point_potential);

#endif // MV

                    if (point_potential > 10.0f)
                    {
                        begin_point_time = actualTime;
                        aux_stim_velocity_flag = true;
                    }
                }
                else
                {
                    point_potential = Vm[0][end_point_index];

#ifdef MV

                    point_potential = rescaleVm(point_potential);

#endif // MV

                    if (point_potential > 10.0f)
                    {
                        end_point_time = actualTime;
                        stim_velocity = (end_point - begin_point) / (end_point_time - begin_point_time); // cm/ms
                        stim_velocity = stim_velocity * 1000.0f;                                         // cm/s
                        stim_velocity_measured = true;
                        INFOMSG("Stim velocity (measured from %.2f to %.2f cm) is %.4g cm/s\n", begin_point, end_point, stim_velocity);
                    }
                }
            }

            finishMeasureVelocityTime = omp_get_wtime();
            elapsedMeasureVelocityTime += finishMeasureVelocityTime - startMeasureVelocityTime;
        }

#endif // MEASURE_VELOCITY

// Update time step counter
#pragma omp single
        timeStepCounter++;
    }

    finishExecutionTime = omp_get_wtime();
    elapsedExecutionTime += finishExecutionTime - startExecutionTime;

#ifdef SAVE_FRAMES

    fprintf(fpFrames, "%lf\n", actualTime);

    saveFrame(fpFrames, Vm, Nx, Ny);

    SUCCESSMSG("Frame at time %.2f ms saved to %s\n", actualTime, framesPath);
    fclose(fpFrames);

#endif // SAVE_FRAMES

    printf("Simulation done!\n");
    printf("\n");

    // Write infos to file
    srand(time(NULL));
    int random_number = (rand() % 500) + 1;
    char infosFilePath[MAX_STRING_SIZE];
    snprintf(infosFilePath, MAX_STRING_SIZE, "%s/infos%d.txt", pathToSaveData, random_number);
    FILE *fpInfos = fopen(infosFilePath, "w");
    fprintf(fpInfos, "EXECUTION TYPE = %s\n", exec_mode);
    fprintf(fpInfos, "PRECISION = %s\n", REAL_TYPE);
    fprintf(fpInfos, "PROBLEM = %s\n", PROBLEM);
    fprintf(fpInfos, "CELL MODEL = %s\n", cell_model);
    fprintf(fpInfos, "METHOD = %s\n", method);
    fprintf(fpInfos, "\n");

    fprintf(fpInfos, "DOMAIN LENGTH IN X = %.4g cm\n", Lx);
    fprintf(fpInfos, "DOMAIN LENGTH IN Y = %.4g cm\n", Ly);
    fprintf(fpInfos, "TOTAL TIME = %.4g ms\n", totalTime);
    fprintf(fpInfos, "\n");

    fprintf(fpInfos, "delta_t = %.5g ms (%d time steps)\n", delta_t, M);
    fprintf(fpInfos, "delta_x = %.5g cm (%d um) (%d space steps in x)\n", delta_x, CM_TO_UM(delta_x), Nx);
    fprintf(fpInfos, "delta_y = %.5g cm (%d um) (%d space steps in y)\n", delta_y, CM_TO_UM(delta_y), Ny);
    fprintf(fpInfos, "TOTAL POINTS IN DOMAIN = %d\n", Nx * Ny);
    fprintf(fpInfos, "DIFFUSION COEFFICIENT = %.8g\n", diff_coeff);
    fprintf(fpInfos, "NUMBER OF STIMULI = %d\n", numberOfStimuli);
    for (int i = 0; i < numberOfStimuli; i++)
    {
        fprintf(fpInfos, "STIMULUS %d: START TIME = %.5g ms\n", i + 1, stimuli[i].begin_time);
    }

    fprintf(fpInfos, "\n");
    fprintf(fpInfos, "NUMBER OF THREADS = %d\n", NUMTHREADS);

#ifdef MEASURE_VELOCITY

    fprintf(fpInfos, "\n");
    fprintf(fpInfos, "STIMULUS S1 VELOCITY = %.4g cm/s\n", stim_velocity);

#endif // MEASURE_VELOCITY

    fprintf(fpInfos, "\n");
    fprintf(fpInfos, "TIME OF THE FIRST PART = %.5g s\n", elapsedTime1stPart);
    fprintf(fpInfos, "TIME OF THE SECOND PART = %.5g s\n", elapsedTime2ndPart);

#if defined(SSIADI) || defined(OSADI)

    fprintf(fpInfos, "TIME TO SOLVE THE 1st LINEAR SYSTEM = %.5g s\n", elapsedTime1stLS);
    fprintf(fpInfos, "TIME TO SOLVE THE 2nd LINEAR SYSTEM = %.5g s\n", elapsedTime2ndLS);

#endif // SSIADI || OSADI

#ifdef MEASURE_VELOCITY

    fprintf(fpInfos, "TIME TO MEASURE VELOCITY = %.5g s\n", elapsedMeasureVelocityTime);

#endif // MEASURE_VELOCITY

#ifdef SAVE_FRAMES

    fprintf(fpInfos, "TIME TO SAVE FRAMES = %.5g s\n", elapsedSaveFramesTime);

#endif // SAVE_FRAMES

    fprintf(fpInfos, "\n");
    fprintf(fpInfos, "SIMULATION TOTAL EXECUTION TIME = %.5g s\n", elapsedExecutionTime);

    fprintf(fpInfos, "\n");
    fprintf(fpInfos, "OUTPUT DIRECTORY = %s\n", pathToSaveData);
    fclose(fpInfos);

    INFOMSG("Simulation total execution time = %.5g s\n", elapsedExecutionTime);
    SUCCESSMSG("Simulation infos saved to %s\n", infosFilePath);

#ifdef SAVE_LAST_FRAME

    // Save last frame
    char lastFrameFilePath[MAX_STRING_SIZE];
    snprintf(lastFrameFilePath, MAX_STRING_SIZE, "%s/lastframe.txt", pathToSaveData);
    FILE *fpLast = fopen(lastFrameFilePath, "w");

    saveFrame(fpLast, Vm, Nx, Ny);

    SUCCESSMSG("Last frame saved to %s\n", lastFrameFilePath);
    fclose(fpLast);

#endif // SAVE_LAST_FRAME

#ifdef SAVE_LAST_STATE

#ifdef AFHN

    char lastFrameFilePathVm[MAX_STRING_SIZE], lastFrameFilePathW[MAX_STRING_SIZE];
    snprintf(lastFrameFilePathVm, MAX_STRING_SIZE, "%s/lastframeVm.txt", pathToSaveData);
    snprintf(lastFrameFilePathW, MAX_STRING_SIZE, "%s/lastframeW.txt", pathToSaveData);

    FILE *fpLastVm = fopen(lastFrameFilePathVm, "w");
    FILE *fpLastW = fopen(lastFrameFilePathW, "w");

    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            fprintf(fpLastVm, "%e ", Vm[i][j]);
            fprintf(fpLastW, "%e ", W[i][j]);
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
    snprintf(lastFrameFilePathVm, MAX_STRING_SIZE, "%s/lastframeVm.txt", pathToSaveData);
    snprintf(lastFrameFilePathX_r1, MAX_STRING_SIZE, "%s/lastframeX_r1.txt", pathToSaveData);
    snprintf(lastFrameFilePathX_r2, MAX_STRING_SIZE, "%s/lastframeX_r2.txt", pathToSaveData);
    snprintf(lastFrameFilePathX_s, MAX_STRING_SIZE, "%s/lastframeX_s.txt", pathToSaveData);
    snprintf(lastFrameFilePathm, MAX_STRING_SIZE, "%s/lastframem.txt", pathToSaveData);
    snprintf(lastFrameFilePathh, MAX_STRING_SIZE, "%s/lastframeh.txt", pathToSaveData);
    snprintf(lastFrameFilePathj, MAX_STRING_SIZE, "%s/lastframej.txt", pathToSaveData);
    snprintf(lastFrameFilePathd, MAX_STRING_SIZE, "%s/lastframed.txt", pathToSaveData);
    snprintf(lastFrameFilePathf, MAX_STRING_SIZE, "%s/lastframef.txt", pathToSaveData);
    snprintf(lastFrameFilePathf2, MAX_STRING_SIZE, "%s/lastframef2.txt", pathToSaveData);
    snprintf(lastFrameFilePathfCaSS, MAX_STRING_SIZE, "%s/lastframefCaSS.txt", pathToSaveData);
    snprintf(lastFrameFilePaths, MAX_STRING_SIZE, "%s/lastframes.txt", pathToSaveData);
    snprintf(lastFrameFilePathr, MAX_STRING_SIZE, "%s/lastframer.txt", pathToSaveData);
    snprintf(lastFrameFilePathR_prime, MAX_STRING_SIZE, "%s/lastframeR_prime.txt", pathToSaveData);
    snprintf(lastFrameFilePathCa_i, MAX_STRING_SIZE, "%s/lastframeCa_i.txt", pathToSaveData);
    snprintf(lastFrameFilePathCa_SR, MAX_STRING_SIZE, "%s/lastframeCa_SR.txt", pathToSaveData);
    snprintf(lastFrameFilePathCa_SS, MAX_STRING_SIZE, "%s/lastframeCa_SS.txt", pathToSaveData);
    snprintf(lastFrameFilePathNa_i, MAX_STRING_SIZE, "%s/lastframeNa_i.txt", pathToSaveData);
    snprintf(lastFrameFilePathK_i, MAX_STRING_SIZE, "%s/lastframeK_i.txt", pathToSaveData);

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
            fprintf(fpLastVm, "%e ", Vm[i][j]);
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

    // Save last state
    char lastStateFilePathVm[MAX_STRING_SIZE], lastStateFilePathv[MAX_STRING_SIZE],
        lastStateFilePathw[MAX_STRING_SIZE], lastStateFilePaths[MAX_STRING_SIZE];
    snprintf(lastStateFilePathVm, MAX_STRING_SIZE, "%s/laststateVm.txt", pathToSaveData);
    snprintf(lastStateFilePathv, MAX_STRING_SIZE, "%s/laststatev.txt", pathToSaveData);
    snprintf(lastStateFilePathw, MAX_STRING_SIZE, "%s/laststatew.txt", pathToSaveData);
    snprintf(lastStateFilePaths, MAX_STRING_SIZE, "%s/laststates.txt", pathToSaveData);

    FILE *fpLastStateVm = fopen(lastStateFilePathVm, "w");
    FILE *fpLastStatev = fopen(lastStateFilePathv, "w");
    FILE *fpLastStatew = fopen(lastStateFilePathw, "w");
    FILE *fpLastStates = fopen(lastStateFilePaths, "w");

    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            fprintf(fpLastStateVm, "%e ", Vm[i][j]);
            fprintf(fpLastStatev, "%e ", v[i][j]);
            fprintf(fpLastStatew, "%e ", w[i][j]);
            fprintf(fpLastStates, "%e ", s[i][j]);
        }
        fprintf(fpLastStateVm, "\n");
        fprintf(fpLastStatev, "\n");
        fprintf(fpLastStatew, "\n");
        fprintf(fpLastStates, "\n");
    }

    SUCCESSMSG("Last Vm state saved to %s\n", lastStateFilePathVm);
    SUCCESSMSG("Last v state saved to %s\n", lastStateFilePathv);
    SUCCESSMSG("Last w state saved to %s\n", lastStateFilePathw);
    SUCCESSMSG("Last s state saved to %s\n", lastStateFilePaths);

    fclose(fpLastStateVm);
    fclose(fpLastStatev);
    fclose(fpLastStatew);
    fclose(fpLastStates);

#endif // MV

#endif // SAVE_LAST_STATE

    // Free memory
    free(time_array);
    free(pathToSaveData);

    for (int i = 0; i < Ny; i++)
    {
        free(Vm[i]);
        free(partRHS[i]);

#if defined(SSIADI)

        free(RHS[i]);

#endif // SSIADI

#ifdef AFHN

        free(W[i]);

#endif // AFHN

#ifdef TT2

        free(X_r1[i]);
        free(X_r2[i]);
        free(X_s[i]);
        free(m[i]);
        free(h[i]);
        free(j[i]);
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

#ifdef MV

        free(v[i]);
        free(w[i]);
        free(s[i]);

#endif // MV
    }

#if defined(SSIADI) || defined(OSADI)

#if defined(SSIADI)

    free(RHS);

#endif // SSIADI

    free(la_y);
    free(lb_y);
    free(lc_y);

#endif // SSIADI || OSADI

    free(Vm);
    free(partRHS);

#if defined(SSIADI) || defined(OSADI)

    free(la_x);
    free(lb_x);
    free(lc_x);

#endif // SSIADI || OSADI

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
    free(R_prime);
    free(Ca_i);
    free(Ca_SR);
    free(Ca_SS);
    free(Na_i);
    free(K_i);

#endif // TT2

#ifdef MV

    free(v);
    free(w);
    free(s);

#endif // MV

    return;
}

#endif // USE_OPENMP