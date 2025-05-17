#include "afhn.h"

#ifndef USE_CUDA

// Model parameters - Based on Gerardo_Giorda 2007
// const real sigma = 1.2e-3f; // omega^-1 * cm^-1
// const real chi = 1.0e3f;    // cm^-1
// const real Cm = 1.0e-3f;    // mF * cm^-2

const real G = 1.5f;      // omega^-1 * cm^-2
const real eta1 = 4.4f;   // omega^-1 * cm^-1
const real eta2 = 0.012f; // dimensionless
const real eta3 = 1.0f;   // dimensionless
const real vth = 13.0f;   // mV
const real vp = 100.0f;   // mV

#endif // USE_CUDA

// RHS functions for the AFHN model
static inline real dVmdt(const real Vm, const real W)
{
    return (G * Vm * (1.0f - (Vm / vth)) * (1.0f - (Vm / vp))) + (eta1 * Vm * W);
}

static inline real dWdt(const real Vm, const real W)
{
    return (eta2 * ((Vm / vp) - (eta3 * W)));
}

// Variables for time measurement
static real startTime = 0.0f;
static real startExecutionTime = 0.0f;
static real elapsedExecutionTime = 0.0f;
static real elapsedTime1stPart = 0.0f;
static real elapsedTime2ndPart = 0.0f;
static real elapsedSaveFramesTime = 0.0f;
static real elapsedMeasureVelocityTime = 0.0f;

// Function to run the simulation using the SSIADI method
void run_SSIADI_AFHN(const SimulationConfig *config, Measurement *measurement, const real *time_array, real *restrict Vm, real *restrict W, real *restrict partRHS)
{
    // Unpack configuration parameters
    const int M = config->M;
    const int Nx = config->Nx;
    const int Ny = config->Ny;
    const real delta_t = config->dt;
    const real delta_x = config->dx;
    const real delta_y = config->dy;
    const real sigma = config->sigma;
    const int numberOfStimuli = config->stimulus_count;
    const Stimulus *stimuli = config->stimuli;

    const bool saveFrames = config->save_frames;
    const int frameSaveRate = config->frame_save_rate;
    const char *pathToSaveData = config->output_dir;
    const char *saveFunctionName = config->save_function_name;
    const char *file_extension = config->file_extension;
    const save_function_t save_function = config->save_function;
    const bool measureVelocity = config->measure_velocity;

    // Extra variables for time measurement
    real startLSTime = 0.0f;
    real elapsedTime1stLS = 0.0f;
    real elapsedTime2ndLS = 0.0f;

    // Measure velocity variables
    real stim_velocity, t0, t1;
    real x0 = config->Lx / 3.0f;
    real x1 = 2.0f * x0;
    int idx_x0 = round(x0 / delta_x) + 1;
    int idx_x1 = round(x1 / delta_x) + 1;
    bool aux_stim_velocity_flag = false;
    bool stim_velocity_measured = false;

    // Auxiliary variables for the loops
    int timeStepCounter = 0;
    real actualTime = 0.0f;
    int i, j, stim_idx;
    int idx, idx_left, idx_right, idx_top, idx_bottom;

    // Auxiliary variables for the operations
    real diff_term, stim, actualVm, actualW;
    real Vmtilde, Wtilde;
    real *RHS = (real *)malloc(Nx * Ny * sizeof(real));

    // Calculate coefficients for the ADI method
    const real phi_x = delta_t / (delta_x * delta_x);
    const real phi_y = delta_t / (delta_y * delta_y);
    const real diff_coeff = sigma / (AFHN_Cm * AFHN_CHI);
    const real tau = 0.5f; // Used for the explicit RHS of ADI

    // Auxiliary arrays for Thomas algorithm
    real *la_x = (real *)malloc(Nx * sizeof(real)); // subdiagonal
    real *lb_x = (real *)malloc(Nx * sizeof(real)); // diagonal
    real *lc_x = (real *)malloc(Nx * sizeof(real)); // superdiagonal
    real *c_prime_x = (real *)malloc(Nx * sizeof(real));
    real *d_prime_x = (real *)malloc(Nx * sizeof(real));
    real *LS_b_x = (real *)malloc(Nx * sizeof(real));
    real *result_x = (real *)malloc(Nx * sizeof(real));
    real *la_y = (real *)malloc(Ny * sizeof(real)); // subdiagonal
    real *lb_y = (real *)malloc(Ny * sizeof(real)); // diagonal
    real *lc_y = (real *)malloc(Ny * sizeof(real)); // superdiagonal
    real *c_prime_y = (real *)malloc(Ny * sizeof(real));
    real *d_prime_y = (real *)malloc(Ny * sizeof(real));
    real *LS_b_y = (real *)malloc(Ny * sizeof(real));
    real *result_y = (real *)malloc(Ny * sizeof(real));

    populateDiagonalThomasAlgorithm(la_x, lb_x, lc_x, Nx, 0.5f * phi_x * diff_coeff);
    populateDiagonalThomasAlgorithm(la_y, lb_y, lc_y, Ny, 0.5f * phi_y * diff_coeff);

    printf("\n");
    printf("Starting simulation...\n");
    startExecutionTime = omp_get_wtime();

    // Main time loop
    while (timeStepCounter < M)
    {
        // Get time step
        actualTime = time_array[timeStepCounter];

        // Start measuring time of 1st part
        startTime = omp_get_wtime();

        // ================================================!
        //  Calculate Approxs. and Update ODEs             !
        // ================================================!
        diff_term = 0.0f;
        for (i = 0; i < Ny; i++)
        {
            for (j = 0; j < Nx; j++)
            {
                idx = i * Nx + j;

                // Calculate the explicit part of the RHS, including the diffusion term in both directions
                actualVm = Vm[idx];
                actualW = W[idx];

                // Stimulation
                stim = get_stimulus_value(actualTime, i, j, stimuli, numberOfStimuli);

                // Calculate aproximation with RK2 -> Vmn+1/2 = Vmn + 0.5*diffusion + 0.5*dt*R(Vmn, Wn)
                diff_term = get_diffusion_term(Vm, i, j, Nx, Ny, diff_coeff, phi_x, phi_y);
                Vmtilde = actualVm + 0.5f * diff_term + (0.5f * delta_t * (stim - dVmdt(actualVm, actualW)));

                // Calculate approximation for state variables and prepare part of the RHS of the following linear systems
                Wtilde = actualW + (0.5f * delta_t * dWdt(actualVm, actualW));
                partRHS[idx] = delta_t * (stim - dVmdt(Vmtilde, Wtilde));

                // Update state variables
                W[idx] = actualW + delta_t * dWdt(Vmtilde, Wtilde); // with RK2 -> Wn+1 = Wn + dt*R(Vm*, W*)
            }
        }

        // End of the 1st part of the time step
        elapsedTime1stPart += omp_get_wtime() - startTime;

        // Start measuring time of 2nd part
        startTime = omp_get_wtime();

        // ================================================!
        //  Calculate Vm at n+1/2 -> Result goes to RHS    !
        //  diffusion implicit in y and explicit in x      !
        // ================================================!
        for (j = 0; j < Nx; j++)
        {
            // Calculate the RHS of the linear system with the explicit diffusion term along x
            for (i = 0; i < Ny; i++)
            {
                idx = i * Nx + j;
                idx_left = i * Nx + lim(j - 1, Nx);
                idx_right = i * Nx + lim(j + 1, Nx);
                actualVm = Vm[idx];
                diff_term = diff_coeff * tau * phi_x * (Vm[idx_left] - 2.0f * actualVm + Vm[idx_right]);
                LS_b_y[i] = actualVm + diff_term + 0.5f * partRHS[idx];
            }

            // Start measuring time of 1st LS
            startLSTime = omp_get_wtime();

            // Solve the linear system
            tridiagonalSystemSolver(la_y, lb_y, lc_y, c_prime_y, d_prime_y, Ny, LS_b_y, result_y);

            // Update with the result
            for (i = 0; i < Ny; i++)
            {
                idx = i * Nx + j;
                RHS[idx] = result_y[i];
            }

            // End measuring time of 1st LS
            elapsedTime1stLS += omp_get_wtime() - startLSTime;
        }

        // ================================================!
        //  Calculate Vm at n+1 -> Result goes to Vm       !
        //  diffusion implicit in x and explicit in y      !
        // ================================================!
        for (i = 0; i < Ny; i++)
        {
            // Calculate the RHS of the linear system with the explicit diffusion term along y
            for (j = 0; j < Nx; j++)
            {
                idx = i * Nx + j;
                idx_top = lim(i + 1, Ny) * Nx + j;
                idx_bottom = lim(i - 1, Ny) * Nx + j;
                actualVm = RHS[idx];
                diff_term = diff_coeff * tau * phi_y * (RHS[idx_bottom] - 2.0f * actualVm + RHS[idx_top]);
                LS_b_x[j] = actualVm + diff_term + 0.5f * partRHS[idx];
            }

            // Start measuring time of 2nd LS
            startLSTime = omp_get_wtime();

            // Solve the linear system
            tridiagonalSystemSolver(la_x, lb_x, lc_x, c_prime_x, d_prime_x, Nx, LS_b_x, result_x);

            // Update with the result
            for (j = 0; j < Nx; j++)
            {
                idx = i * Nx + j;
                Vm[idx] = result_x[j];
            }

            // End measuring time of 2nd LS
            elapsedTime2ndLS += omp_get_wtime() - startLSTime;
        }

        // End of the 2nd part of the time step
        elapsedTime2ndPart += omp_get_wtime() - startTime;

        // Save frame if needed
        startTime = omp_get_wtime();
        handle_frame_saving(pathToSaveData, file_extension, save_function, timeStepCounter, frameSaveRate, saveFrames, Vm, Nx, Ny, delta_x, delta_y, actualTime);
        elapsedSaveFramesTime += omp_get_wtime() - startTime;

        // Measure velocity if needed
        startTime = omp_get_wtime();
        handle_velocity_measurement(Vm, idx_x0, idx_x1, &t0, &t1, &aux_stim_velocity_flag, &stim_velocity_measured, actualTime, x0, x1, &stim_velocity);
        elapsedMeasureVelocityTime += omp_get_wtime() - startTime;

        // Update time step counter
        timeStepCounter++;
    }

    elapsedExecutionTime += omp_get_wtime() - startExecutionTime;

    printf("Simulation done!\n");
    printf("\n");

    // Update measurement structure
    measurement->elapsedExecutionTime = elapsedExecutionTime;
    measurement->elapsedTime1stPart = elapsedTime1stPart;
    measurement->elapsedTime2ndPart = elapsedTime2ndPart;
    measurement->elapsedTime1stLS = elapsedTime1stLS;
    measurement->elapsedTime2ndLS = elapsedTime2ndLS;
    measurement->elapsedSaveFramesTime = elapsedSaveFramesTime;
    measurement->elapsedMeasureVelocityTime = elapsedMeasureVelocityTime;
    measurement->stimVelocity = stim_velocity;

    // Free allocated memory
    free(RHS);
    free(la_x);
    free(lb_x);
    free(lc_x);
    free(c_prime_x);
    free(d_prime_x);
    free(LS_b_x);
    free(result_x);
    free(la_y);
    free(lb_y);
    free(lc_y);
    free(c_prime_y);
    free(d_prime_y);
    free(LS_b_y);
    free(result_y);
}

// Function to run the simulation using the OSADI method
void run_OSADI_AFHN(const SimulationConfig *config, Measurement *measurement, const real *time_array, real *restrict Vm, real *restrict W, real *restrict partRHS)
{
    // Unpack configuration parameters
    const int M = config->M;
    const int Nx = config->Nx;
    const int Ny = config->Ny;
    const real delta_t = config->dt;
    const real delta_x = config->dx;
    const real delta_y = config->dy;
    const real sigma = config->sigma;
    const int numberOfStimuli = config->stimulus_count;
    const Stimulus *stimuli = config->stimuli;

    const bool saveFrames = config->save_frames;
    const int frameSaveRate = config->frame_save_rate;
    const char *pathToSaveData = config->output_dir;
    const char *saveFunctionName = config->save_function_name;
    const char *file_extension = config->file_extension;
    const save_function_t save_function = config->save_function;
    const bool measureVelocity = config->measure_velocity;

    // Variables for time measurement
    real startLSTime = 0.0f;
    real elapsedTime1stLS = 0.0f;
    real elapsedTime2ndLS = 0.0f;

    // Measure velocity variables
    real stim_velocity, t0, t1;
    real x0 = config->Lx / 3.0f;
    real x1 = 2.0f * x0;
    int idx_x0 = round(x0 / delta_x) + 1;
    int idx_x1 = round(x1 / delta_x) + 1;
    bool aux_stim_velocity_flag = false;
    bool stim_velocity_measured = false;

    // Auxiliary variables for the loops
    int timeStepCounter = 0;
    real actualTime = 0.0f;
    int i, j, stim_idx;
    int idx, idx_left, idx_right, idx_top, idx_bottom;

    // Auxiliary variables for the operations
    real stim, diff_term, actualVm, actualW;

    // Calculate coefficients for the ADI method
    const real phi_x = delta_t / (delta_x * delta_x);
    const real phi_y = delta_t / (delta_y * delta_y);
    const real diff_coeff = sigma / (AFHN_Cm * AFHN_CHI);

    // Auxiliary arrays for Thomas algorithm
    real *la_x = (real *)malloc(Nx * sizeof(real)); // subdiagonal
    real *lb_x = (real *)malloc(Nx * sizeof(real)); // diagonal
    real *lc_x = (real *)malloc(Nx * sizeof(real)); // superdiagonal
    real *c_prime_x = (real *)malloc(Nx * sizeof(real));
    real *d_prime_x = (real *)malloc(Nx * sizeof(real));
    real *LS_b_x = (real *)malloc(Nx * sizeof(real));
    real *result_x = (real *)malloc(Nx * sizeof(real));
    real *la_y = (real *)malloc(Ny * sizeof(real)); // subdiagonal
    real *lb_y = (real *)malloc(Ny * sizeof(real)); // diagonal
    real *lc_y = (real *)malloc(Ny * sizeof(real)); // superdiagonal
    real *c_prime_y = (real *)malloc(Ny * sizeof(real));
    real *d_prime_y = (real *)malloc(Ny * sizeof(real));
    real *LS_b_y = (real *)malloc(Ny * sizeof(real));
    real *result_y = (real *)malloc(Ny * sizeof(real));

    populateDiagonalThomasAlgorithm(la_x, lb_x, lc_x, Nx, phi_x * diff_coeff);
    populateDiagonalThomasAlgorithm(la_y, lb_y, lc_y, Ny, phi_y * diff_coeff);

    printf("\n");
    printf("Starting simulation...\n");
    startExecutionTime = omp_get_wtime();

    // Main time loop
    while (timeStepCounter < M)
    {
        // Get time step
        actualTime = time_array[timeStepCounter];

        // Start measuring time of 1st part
        startTime = omp_get_wtime();

        // ================================================!
        //  Calculate Approxs. and Update ODEs             !
        // ================================================!
        diff_term = 0.0f;
        for (i = 0; i < Ny; i++)
        {
            for (j = 0; j < Nx; j++)
            {
                idx = i * Nx + j;

                // Calculate the explicit part of the RHS, including the diffusion term in both directions
                actualVm = Vm[idx];
                actualW = W[idx];

                // Stimulation
                stim = get_stimulus_value(actualTime, i, j, stimuli, numberOfStimuli);

                // Calculate part of the RHS of the following linear systems with Forward Euler
                partRHS[idx] = delta_t * (stim - dVmdt(actualVm, actualW));

                // Update state variables
                W[idx] = actualW + delta_t * dWdt(actualVm, actualW); // with Forward Euler -> Wn+1 = Wn + dt*R(Vmn, Wn)
            }
        }

        // End of the 1st part of the time step
        elapsedTime1stPart += omp_get_wtime() - startTime;

        // Start measuring time of 2nd part
        startTime = omp_get_wtime();

        // ================================================!
        //  Calculate Vm at n+1/2 -> Result goes to Vm       !
        // ================================================!
        for (j = 0; j < Nx; j++)
        {
            // Calculate the RHS of the linear system
            for (i = 0; i < Ny; i++)
            {
                idx = i * Nx + j;
                LS_b_y[i] = Vm[idx] + 0.5f * partRHS[idx];
            }

            // Start measuring time of 1st LS
            startLSTime = omp_get_wtime();

            // Solve the linear system
            tridiagonalSystemSolver(la_y, lb_y, lc_y, c_prime_y, d_prime_y, Ny, LS_b_y, result_y);

            // Update with the result
            for (i = 0; i < Ny; i++)
            {
                idx = i * Nx + j;
                Vm[idx] = result_y[i];
            }

            // End measuring time of 1st LS
            elapsedTime1stLS += omp_get_wtime() - startLSTime;
        }

        // ================================================!
        //  Calculate Vm at n+1 -> Result goes to Vm         !
        // ================================================!
        for (i = 0; i < Ny; i++)
        {
            // Calculate the RHS of the linear system
            for (j = 0; j < Nx; j++)
            {
                idx = i * Nx + j;
                LS_b_x[j] = Vm[idx] + 0.5f * partRHS[idx];
            }

            // Start measuring time of 2nd LS
            startLSTime = omp_get_wtime();

            // Solve the linear system
            tridiagonalSystemSolver(la_x, lb_x, lc_x, c_prime_x, d_prime_x, Nx, LS_b_x, result_x);

            // Update with the result
            for (j = 0; j < Nx; j++)
            {
                idx = i * Nx + j;
                Vm[idx] = result_x[j];
            }

            // End measuring time of 2nd LS
            elapsedTime2ndLS += omp_get_wtime() - startLSTime;
        }

        // End of the 2nd part of the time step
        elapsedTime2ndPart += omp_get_wtime() - startTime;

        // Save frame if needed
        startTime = omp_get_wtime();
        handle_frame_saving(pathToSaveData, file_extension, save_function, timeStepCounter, frameSaveRate, saveFrames, Vm, Nx, Ny, delta_x, delta_y, actualTime);
        elapsedSaveFramesTime += omp_get_wtime() - startTime;

        // Measure velocity if needed
        startTime = omp_get_wtime();
        handle_velocity_measurement(Vm, idx_x0, idx_x1, &t0, &t1, &aux_stim_velocity_flag, &stim_velocity_measured, actualTime, x0, x1, &stim_velocity);
        elapsedMeasureVelocityTime += omp_get_wtime() - startTime;

        // Update time step counter
        timeStepCounter++;
    }

    elapsedExecutionTime += omp_get_wtime() - startExecutionTime;

    printf("Simulation done!\n");
    printf("\n");

    // Update measurement structure
    measurement->elapsedExecutionTime = elapsedExecutionTime;
    measurement->elapsedTime1stPart = elapsedTime1stPart;
    measurement->elapsedTime2ndPart = elapsedTime2ndPart;
    measurement->elapsedTime1stLS = elapsedTime1stLS;
    measurement->elapsedTime2ndLS = elapsedTime2ndLS;
    measurement->elapsedSaveFramesTime = elapsedSaveFramesTime;
    measurement->elapsedMeasureVelocityTime = elapsedMeasureVelocityTime;
    measurement->stimVelocity = stim_velocity;

    // Free allocated memory
    free(la_x);
    free(lb_x);
    free(lc_x);
    free(c_prime_x);
    free(d_prime_x);
    free(LS_b_x);
    free(result_x);
    free(la_y);
    free(lb_y);
    free(lc_y);
    free(c_prime_y);
    free(d_prime_y);
    free(LS_b_y);
    free(result_y);
}

// Function to run the simulation using the FE method
void run_FE_AFHN(const SimulationConfig *config, Measurement *measurement, const real *time_array, real *restrict Vm, real *restrict W, real *restrict partRHS)
{
    // Unpack configuration parameters
    const int M = config->M;
    const int Nx = config->Nx;
    const int Ny = config->Ny;
    const real delta_t = config->dt;
    const real delta_x = config->dx;
    const real delta_y = config->dy;
    const real sigma = config->sigma;
    const int numberOfStimuli = config->stimulus_count;
    const Stimulus *stimuli = config->stimuli;

    const bool saveFrames = config->save_frames;
    const int frameSaveRate = config->frame_save_rate;
    const char *pathToSaveData = config->output_dir;
    const char *saveFunctionName = config->save_function_name;
    const char *file_extension = config->file_extension;
    const save_function_t save_function = config->save_function;
    const bool measureVelocity = config->measure_velocity;

    // Measure velocity variables
    real stim_velocity, t0, t1;
    real x0 = config->Lx / 3.0f;
    real x1 = 2.0f * x0;
    int idx_x0 = round(x0 / delta_x) + 1;
    int idx_x1 = round(x1 / delta_x) + 1;
    bool aux_stim_velocity_flag = false;
    bool stim_velocity_measured = false;

    // Auxiliary variables for the loops
    int timeStepCounter = 0;
    real actualTime = 0.0f;
    int i, j, stim_idx;
    int idx, idx_left, idx_right, idx_top, idx_bottom;

    // Auxiliary variables for the operations
    real stim, diff_term, actualVm, actualW;

    // Calculate coefficients for the ADI method
    const real phi_x = delta_t / (delta_x * delta_x);
    const real phi_y = delta_t / (delta_y * delta_y);
    const real diff_coeff = sigma / (AFHN_Cm * AFHN_CHI);

    printf("\n");
    printf("Starting simulation...\n");
    startExecutionTime = omp_get_wtime();

    // Main time loop
    while (timeStepCounter < M)
    {
        // Get time step
        actualTime = time_array[timeStepCounter];

        // Start measuring time of 1st part
        startTime = omp_get_wtime();

        // ================================================!
        //  Calculate Approxs. and Update ODEs             !
        // ================================================!
        diff_term = 0.0f;
        for (i = 0; i < Ny; i++)
        {
            for (j = 0; j < Nx; j++)
            {
                idx = i * Nx + j;

                // Calculate the explicit part of the RHS, including the diffusion term in both directions
                actualVm = Vm[idx];
                actualW = W[idx];

                // Stimulation
                stim = get_stimulus_value(actualTime, i, j, stimuli, numberOfStimuli);

                // Update variables explicitly
                diff_term = get_diffusion_term(Vm, i, j, Nx, Ny, diff_coeff, phi_x, phi_y);
                partRHS[idx] = actualVm + diff_term + delta_t * (stim - dVmdt(actualVm, actualW));

                W[idx] = actualW + delta_t * dWdt(actualVm, actualW);
            }
        }

        // End of the 1st part of the time step
        elapsedTime1stPart += omp_get_wtime() - startTime;

        // Start measuring time of 2nd part
        startTime = omp_get_wtime();

        // ==================!
        //  Update Vm        !
        // ==================!
        for (idx = 0; idx < Nx * Ny; idx++)
            Vm[idx] = partRHS[idx];

        // End of the 2nd part of the time step
        elapsedTime2ndPart += omp_get_wtime() - startTime;

        // Save frame if needed
        startTime = omp_get_wtime();
        handle_frame_saving(pathToSaveData, file_extension, save_function, timeStepCounter, frameSaveRate, saveFrames, Vm, Nx, Ny, delta_x, delta_y, actualTime);
        elapsedSaveFramesTime += omp_get_wtime() - startTime;

        // Measure velocity if needed
        startTime = omp_get_wtime();
        handle_velocity_measurement(Vm, idx_x0, idx_x1, &t0, &t1, &aux_stim_velocity_flag, &stim_velocity_measured, actualTime, x0, x1, &stim_velocity);
        elapsedMeasureVelocityTime += omp_get_wtime() - startTime;

        // Update time step counter
        timeStepCounter++;
    }

    elapsedExecutionTime += omp_get_wtime() - startExecutionTime;

    printf("Simulation done!\n");
    printf("\n");

    // Update measurement structure
    measurement->elapsedExecutionTime = elapsedExecutionTime;
    measurement->elapsedTime1stPart = elapsedTime1stPart;
    measurement->elapsedTime2ndPart = elapsedTime2ndPart;
    measurement->elapsedSaveFramesTime = elapsedSaveFramesTime;
    measurement->elapsedMeasureVelocityTime = elapsedMeasureVelocityTime;
    measurement->stimVelocity = stim_velocity;
}

// Map of numerical methods and their corresponding solvers for AFHN model
static const struct
{
    const NumericalMethod method;
    numerical_method_afhn_t solver;
} numerical_method_afhn_map[] = {
    {METHOD_SSIADI, run_SSIADI_AFHN},
    {METHOD_OSADI, run_OSADI_AFHN},
    {METHOD_FE, run_FE_AFHN},
    {METHOD_INVALID, NULL}};

numerical_method_afhn_t get_numerical_method_afhn(const NumericalMethod *method)
{
    if (method == NULL)
    {
        printf("Error: numerical method is NULL\n");
        return NULL;
    }

    for (int i = 0; numerical_method_afhn_map[i].solver != NULL; i++)
        if (numerical_method_afhn_map[i].method == *method)
            return numerical_method_afhn_map[i].solver;
    return NULL;
}

real forcingTerm(real x, real y, real t, real W, real Lx, real Ly, real sigma)
{
    // Calculate coefficients for the ADI method
    const real diff_coeff = sigma / (AFHN_Cm * AFHN_CHI);
    real exactVm = (exp(-t)) * cos(_PI * x / Lx) * cos(_PI * y / Ly);
    real reaction = dVmdt(exactVm, W);
    return (exactVm * (-(AFHN_CHI * AFHN_Cm) + 2.0f * (sigma / (AFHN_CHI * AFHN_Cm)) * _PI * _PI / (Lx * Ly))) + (AFHN_CHI * reaction);
}

void initializeWithInitialConditionAFHN(const real Nx, const real Ny, real *Vm, real *W)
{
    // Initial conditions
    const real Vm_init = 0.0f;
    const real W_init = 0.0f;

    // Initialize Vm and W with the initial condition
    for (int idx = 0; idx < Nx * Ny; idx++)
    {
        Vm[idx] = Vm_init;
        W[idx] = W_init;
    }
}

void solveMonodomainAFHN(const SimulationConfig *config, Measurement *measurement, const real *time_array)
{
    // Allocate aux array for variables
    int total_points = config->Nx * config->Ny;
    real *partRHS = (real *)malloc(total_points * sizeof(real));

    // Allocate and initialize state variables arrays
    real *Vm = (real *)malloc(total_points * sizeof(real));
    real *W = (real *)malloc(total_points * sizeof(real));
    initializeWithInitialConditionAFHN(config->Nx, config->Ny, Vm, W);

    if (config->path_to_restore_state_files != NULL && strlen(config->path_to_restore_state_files) > 0)
    {
        printf("\n");
        printf("Restoring state variables...\n");

        // Load initial conditions from file
        // TODO: Implement file loading logic
    }

    if (config->shift_state)
    {
        printf("\n");
        printf("Shifting state variables...\n");

        // Shift state variables to the left
        real lengthToShift = 0.0f;
        // TODO: Implement shift logic
    }

    // Run the simulation based on the selected method
    numerical_method_afhn_t run_method = get_numerical_method_afhn(&config->method);
    if (run_method == NULL)
    {
        fprintf(stderr, "Error: Invalid numerical method for AFHN model\n");
        free(partRHS);
        free(Vm);
        free(W);
        return;
    }

    // Run the selected method
    run_method(config, measurement, time_array, Vm, W, partRHS);

    // Save last frame
    if (config->save_last_frame)
    {
        startTime = omp_get_wtime();
        static char file_path[MAX_STRING_SIZE];
        snprintf(file_path, MAX_STRING_SIZE, "%s/frames/Vm_%05d.%s", config->output_dir, config->M, config->file_extension);
        config->save_function(file_path, Vm, config->Nx, config->Ny, config->dx, config->dy);
        measurement->elapsedSaveFramesTime += (omp_get_wtime() - startTime);
        SUCCESSMSG("Last frame (time %.2f ms) saved to %s\n", config->M * config->dt, file_path);
    }

    // Save last state
    if (config->save_last_state)
    {
        startTime = omp_get_wtime();
        // TODO: save as binary
        static char file_path[MAX_STRING_SIZE];
        snprintf(file_path, MAX_STRING_SIZE, "%s/state_%05d.dat", config->output_dir, config->M);
        measurement->elapsedSaveStateTime += (omp_get_wtime() - startTime);
        SUCCESSMSG("Last state (time %.2f ms) saved to %s\n", config->M * config->dt, file_path);
    }

    // Free memory
    free(partRHS);
    free(Vm);
    free(W);
}



// Instantiate the AFHN model
const CellModelSolver AFHN_MODEL = {
    .n_state_vars = AFHN_NSV,
    .initialize = initialize_AFHN,
    .compute_dVmdt = compute_dVmdt_AFHN,
    .compute_dSdt = compute_dSdt_AFHN,
    .update_sV = update_sV_AFHN,
};