#include "../numerical_methods.h"
#include "../numerical_methods_helpers.h"

void runOSADI(const SimulationConfig *config, Measurement *measurement, const real *restrict time_array, const CellModelSolver *cell_model_solver, real *Vm, real *sV, real *partRHS)
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
    const int n_state_vars = cell_model_solver->n_state_vars;
    const get_actual_sV_t get_actual_sV = cell_model_solver->get_actual_sV;
    const compute_dVmdt_t compute_dVmdt = cell_model_solver->compute_dVmdt;
    const update_sV_t update_sV = cell_model_solver->update_sV;

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
    real *actualsV = (real *)malloc(cell_model_solver->n_state_vars * sizeof(real));
    if (actualsV == NULL)
    {
        fprintf(stderr, "Error allocating memory for actualsV\n");
        exit(EXIT_FAILURE);
    }
    real *dSdt = (real *)malloc(cell_model_solver->n_state_vars * sizeof(real));
    if (dSdt == NULL)
    {
        fprintf(stderr, "Error allocating memory for dSdt\n");
        exit(EXIT_FAILURE);
    }

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

    // Variables for time measurement
    real startTime = 0.0f;
    real startIterationTime = 0.0f;
    real elapsedExecutionTime = 0.0f;
    real elapsedTime1stPart = 0.0f;
    real elapsedTime2ndPart = 0.0f;
    real elapsedSaveFramesTime = 0.0f;
    real elapsedMeasureVelocityTime = 0.0f;
    real startLSTime = 0.0f;
    real elapsedTime1stLS = 0.0f;
    real elapsedTime2ndLS = 0.0f;

    printf("\n");
    printf("Starting simulation...\n");

    // Main time loop
    while (timeStepCounter < M)
    {
        startIterationTime = omp_get_wtime();

        // Get time step
        actualTime = time_array[timeStepCounter];

        // ================================================!
        //  Calculate Approxs. and Update ODEs             !
        // ================================================!
        startTime = omp_get_wtime();

        diff_term = 0.0f;
        for (i = 0; i < Ny; i++)
        {
            for (j = 0; j < Nx; j++)
            {
                idx = i * Nx + j;

                // Get the actual Vm and sV
                actualVm = Vm[idx];
                get_actual_sV(actualsV, sV, idx);

                // Stimulation
                stim = get_stimulus_value(actualTime, i, j, stimuli, numberOfStimuli);

                // Calculate part of the RHS of the following linear systems with Forward Euler
                partRHS[idx] = delta_t * (stim - compute_dVmdt(actualVm, actualsV));

                // Update state variables
                update_sV(sV, actualsV, actualVm, actualsV, delta_t, idx);
            }
        }

        elapsedTime1stPart += omp_get_wtime() - startTime;

        // ================================================!
        //  Calculate Vm at n+1/2 -> Result goes to Vm     !
        // ================================================!
        startTime = omp_get_wtime();

        for (j = 0; j < Nx; j++)
        {
            // Calculate the RHS of the linear system
            for (i = 0; i < Ny; i++)
            {
                idx = i * Nx + j;
                LS_b_y[i] = Vm[idx] + 0.5f * partRHS[idx];
            }

            // Solve the linear system
            startLSTime = omp_get_wtime();

            tridiagonalSystemSolver(la_y, lb_y, lc_y, c_prime_y, d_prime_y, Ny, LS_b_y, result_y);

            // Update with the result
            for (i = 0; i < Ny; i++)
            {
                idx = i * Nx + j;
                Vm[idx] = result_y[i];
            }

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

            // Solve the linear system
            startLSTime = omp_get_wtime();

            tridiagonalSystemSolver(la_x, lb_x, lc_x, c_prime_x, d_prime_x, Nx, LS_b_x, result_x);

            // Update with the result
            for (j = 0; j < Nx; j++)
            {
                idx = i * Nx + j;
                Vm[idx] = result_x[j];
            }

            elapsedTime2ndLS += omp_get_wtime() - startLSTime;
        }

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

        elapsedExecutionTime += omp_get_wtime() - startIterationTime;
    }

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
    free(actualsV);
    free(dSdt);
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