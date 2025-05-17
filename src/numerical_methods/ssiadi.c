#include "numerical_methods.h"
#include "numerical_methods_helpers.h"

// Variables for time measurement
static real startTime = 0.0f;
static real startExecutionTime = 0.0f;
static real elapsedExecutionTime = 0.0f;
static real elapsedTime1stPart = 0.0f;
static real elapsedTime2ndPart = 0.0f;
static real elapsedSaveFramesTime = 0.0f;
static real elapsedMeasureVelocityTime = 0.0f;

void runSSIADI(const SimulationConfig *config, Measurement *measurement, const real *restrict time_array, const CellModelSolver *cell_model_solver, real *restrict Vm, real *restrict sV, real *restrict partRHS)
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

    // Get the solver functions
    const compute_dVmdt_t compute_dVmdt = cell_model_solver->compute_dVmdt;
    const compute_dSdt_t compute_dSdt = cell_model_solver->compute_dSdt;
    const update_sV_t update_sV = cell_model_solver->update_sV;

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
    real diff_term, stim, actualVm;
    real Vmtilde;
    real *dSdt = (real *)malloc(cell_model_solver->n_state_vars * sizeof(real));
    if (dSdt == NULL)
    {
        fprintf(stderr, "Error allocating memory for dSdt\n");
        exit(EXIT_FAILURE);
    }
    real *sVtilde = (real *)malloc(Nx * Ny * cell_model_solver->n_state_vars * sizeof(real));
    if (sVtilde == NULL)
    {
        fprintf(stderr, "Error allocating memory for sVtilde\n");
        exit(EXIT_FAILURE);
    }
    real *RHS = (real *)malloc(Nx * Ny * sizeof(real));
    if (RHS == NULL)
    {
        fprintf(stderr, "Error allocating memory for RHS\n");
        exit(EXIT_FAILURE);
    }

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

                // Stimulation
                stim = get_stimulus_valueNM(actualTime, i, j, stimuli, numberOfStimuli);

                // Calculate aproximation with RK2 -> Vmn+1/2 = Vmn + 0.5*diffusion + 0.5*dt*R(Vmn, Wn)
                diff_term = get_diffusion_termNM(Vm, i, j, Nx, Ny, diff_coeff, phi_x, phi_y);
                Vmtilde = actualVm + 0.5f * diff_term + (0.5f * delta_t * (stim - compute_dVmdt(actualVm, sV, idx)));

                // Calculate approximation for state variables and prepare part of the RHS of the following linear systems
                compute_dSdt(actualVm, sV, dSdt, idx);
                update_sV(sVtilde, sV, dSdt, 0.5f * delta_t, idx);
                partRHS[idx] = delta_t * (stim - compute_dVmdt(Vmtilde, sVtilde, idx));

                // Update state variables
                compute_dSdt(Vmtilde, sVtilde, dSdt, idx);
                update_sV(sV, sV, dSdt, delta_t, idx);
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
        handle_frame_savingNM(pathToSaveData, file_extension, save_function, timeStepCounter, frameSaveRate, saveFrames, Vm, Nx, Ny, delta_x, delta_y, actualTime);
        elapsedSaveFramesTime += omp_get_wtime() - startTime;

        // Measure velocity if needed
        startTime = omp_get_wtime();
        handle_velocity_measurementNM(Vm, idx_x0, idx_x1, &t0, &t1, &aux_stim_velocity_flag, &stim_velocity_measured, actualTime, x0, x1, &stim_velocity);
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
    free(dSdt);
    free(sVtilde);
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