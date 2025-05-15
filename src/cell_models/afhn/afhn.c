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

real forcingTerm(real x, real y, real t, real W, real Lx, real Ly, real sigma)
{
    // Calculate coefficients for the ADI method
    const real chi = 1.0e3f;    // cm^-1
    const real Cm = 1.0e-3f;    // mF * cm^-2
    const real diff_coeff = sigma / (Cm * chi);
    real exactVm = (exp(-t)) * cos(_PI * x / Lx) * cos(_PI * y / Ly);
    real reaction = (G * exactVm * (1.0f - (exactVm / vth)) * (1.0f - (exactVm / vp))) + (eta1 * exactVm * W);
    return (exactVm * (-(chi * Cm) + 2.0f * (sigma / (chi * Cm)) * _PI * _PI / (Lx * Ly))) + (chi * reaction);
}

void initializeWithInitialConditionAFHN(const real Nx, const real Ny, real *Vm, real *W)
{
    // Initial conditions
    const real Vm_init = 0.0f;
    const real W_init = 0.0f;

    // Initialize Vm and W with the initial condition
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            int index = i * Nx + j;
            Vm[index] = Vm_init;
            W[index] = W_init;
        }
    }
}

void solveMonodomainAFHN(const SimulationConfig *config, Measurement *measurement, const real *time_array)
{
    // Unpack configuration parameters
    const int M = config->M; // Number of time steps
    const int Nx = config->Nx; // Spatial steps in x
    const int Ny = config->Ny; // Spatial steps in y

    // Allocate aux array for variables
    int total_points = Nx * Ny;
    real *partRHS = (real *)malloc(total_points * sizeof(real));

    // Allocate and initialize state variables arrays
    real *Vm = (real *)malloc(total_points * sizeof(real));
    real *W = (real *)malloc(total_points * sizeof(real));
    initializeWithInitialConditionAFHN(Nx, Ny, Vm, W);

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
    if (config->method == METHOD_SSIADI)
    {
        // Run the simulation using the SSIADI method
        run_SSIADI_AFHN(config, measurement, time_array, Vm, W, partRHS);
    }
    else if (config->method == METHOD_OSADI)
    {
        // Run the simulation using the OSADI method
        run_OSADI_AFHN(config, measurement, time_array, Vm, W, partRHS);
    }
    else if (config->method == METHOD_FE)
    {
        // Run the simulation using the Forward Euler method
        run_FE_AFHN(config, measurement, time_array, Vm, W, partRHS);
    }
    else
    {
        fprintf(stderr, "Error: Unsupported method for AFHN model.\n");
        exit(EXIT_FAILURE);
    }

    // Save frames variables
    char file_path[MAX_STRING_SIZE];
    char file_extension[4];
    if (strstr(config->save_function_name, "txt") != NULL)
        snprintf(file_extension, 4, "txt");
    else if (strstr(config->save_function_name, "vtk") != NULL)
        snprintf(file_extension, 4, "vtk");

    // Save last frame
    if (config->save_last_state)
    {
        snprintf(file_path, MAX_STRING_SIZE, "%s/Vm_%05d.%s", config->output_dir, M, file_extension);
        config->save_function(file_path, Vm, Nx, Ny, config->dx, config->dy);
        SUCCESSMSG("Last frame (time %.2f ms) saved to %s\n", M * config->dt, file_path);
    }

    // Save last state
    if (config->save_last_state)
    {
        // TODO: save as binary
        snprintf(file_path, MAX_STRING_SIZE, "%s/state_%05d.dat", config->output_dir, M);
        SUCCESSMSG("Last state (time %.2f ms) saved to %s\n", M * config->dt, file_path);
    }

    // Free memory
    free(partRHS);
    free(Vm);
    free(W);
}

// Function to run the simulation using the SSIADI method
void run_SSIADI_AFHN(const SimulationConfig *config, Measurement *measurement, const real *time_array, real *Vm, real *W, real *partRHS)
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
    const save_function_t save_function = config->save_function;
    const bool measureVelocity = config->measure_velocity;

    // Variables for time measurement
    real startTime = 0.0f;
    real finishTime = 0.0f;
    real startExecutionTime = 0.0f;
    real finishExecutionTime = 0.0f;
    real elapsedExecutionTime = 0.0f;
    real elapsedTime1stPart = 0.0f;
    real elapsedTime2ndPart = 0.0f;
    real startLSTime = 0.0f;
    real finishLSTime = 0.0f;
    real elapsedTime1stLS = 0.0f;
    real elapsedTime2ndLS = 0.0f;
    real startSaveFramesTime = 0.0f;
    real finishSaveFramesTime = 0.0f;
    real elapsedSaveFramesTime = 0.0f;
    real startMeasureVelocityTime = 0.0f;
    real finishMeasureVelocityTime = 0.0f;
    real elapsedMeasureVelocityTime = 0.0f;

    // Auxiliary variables for the loops
    int timeStepCounter = 0;
    real actualTime = 0.0f;
    int i, j, index, si;

    // Save frames variables
    char file_path[MAX_STRING_SIZE];
    char file_extension[4];
    if (strstr(saveFunctionName, "txt") != NULL)
        snprintf(file_extension, 4, "txt");
    else if (strstr(saveFunctionName, "vtk") != NULL)
        snprintf(file_extension, 4, "vtk");

    // Measure velocity variables
    real stim_velocity = 0.0f;
    real point_potential = 0.0f;
    real begin_point = config->Lx / 3.0f;
    real end_point = 2.0f * begin_point;
    int begin_point_index = round(begin_point / delta_x) + 1;
    int end_point_index = round(end_point / delta_x) + 1;
    real begin_point_time = 0.0f;
    real end_point_time = 0.0f;
    bool aux_stim_velocity_flag = false;
    bool stim_velocity_measured = false;

    // Auxiliary variables for the operations
    real *RHS = (real *)malloc(Nx * Ny * sizeof(real));
    real stim = 0.0f;
    real diff_term = 0.0f;
    real actualVm;
    real actualW;
    real RHS_Vm_term;
    real Vmtilde;
    real Wtilde;
    real RHS_Vmtilde_term;

    // Calculate coefficients for the ADI method
    const real phi_x = delta_t / (delta_x * delta_x);
    const real phi_y = delta_t / (delta_y * delta_y);
    const real tau = 0.5f; // Used for calculating the explicit diffusion term on the rhs of the SSIADI method
    const real chi = 1.0e3f;    // cm^-1
    const real Cm = 1.0e-3f;    // mF * cm^-2
    const real diff_coeff = sigma / (Cm * chi);

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
                index = i * Nx + j;

                // Calculate the explicit part of the RHS, including the diffusion term in both directions
                actualVm = Vm[index];
                actualW = W[index];
                RHS_Vm_term = (G * actualVm * (1.0f - (actualVm / vth)) * (1.0f - (actualVm / vp))) + (eta1 * actualVm * actualW);

                // Stimulation
                stim = 0.0f;

#pragma unroll
                for (si = 0; si < numberOfStimuli; si++)
                {
                    if (actualTime >= stimuli[si].begin_time && actualTime <= stimuli[si].begin_time + stimuli[si].duration && j >= stimuli[si].x_discretized.min && j <= stimuli[si].x_discretized.max && i >= stimuli[si].y_discretized.min && i <= stimuli[si].y_discretized.max)
                    {
                        stim = stimuli[si].amplitude;
                        break;
                    }
                }

                // Calculate aproximation with RK2 -> Vmn+1/2 = Vmn + 0.5*diffusion + 0.5*dt*R(Vmn, Wn)
                diff_term = diff_coeff * (phi_x * (Vm[i * Nx + (lim(j - 1, Nx))] - 2.0f * actualVm + Vm[i * Nx + lim(j + 1, Nx)]) + phi_y * (Vm[(lim(i - 1, Ny) * Nx + j)] - 2.0f * actualVm + Vm[(lim(i + 1, Ny) * Nx + j)]));
                Vmtilde = actualVm + 0.5f * diff_term + (0.5f * delta_t * (stim - RHS_Vm_term));

                // Calculate approximation for state variables and prepare part of the RHS of the following linear systems
                Wtilde = actualW + (0.5f * delta_t * (eta2 * ((actualVm / vp) - (eta3 * actualW))));
                RHS_Vmtilde_term = (G * Vmtilde * (1.0f - (Vmtilde / vth)) * (1.0f - (Vmtilde / vp))) + (eta1 * Vmtilde * Wtilde);
                partRHS[index] = delta_t * (stim - RHS_Vmtilde_term);

                // Update state variables
                W[index] = actualW + delta_t * (eta2 * ((Vmtilde / vp) - (eta3 * Wtilde))); // with RK2 -> Wn+1 = Wn + dt*R(Vm*, W*)
            }
        }

        // End of the 1st part of the time step
        finishTime = omp_get_wtime();
        elapsedTime1stPart += finishTime - startTime;

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
                index = i * Nx + j;
                actualVm = Vm[index];
                diff_term = diff_coeff * tau * phi_x * (Vm[i * Nx + (lim(j - 1, Nx))] - 2.0f * actualVm + Vm[i * Nx + lim(j + 1, Nx)]);
                LS_b_y[i] = actualVm + diff_term + 0.5f * partRHS[index];
            }

            // Start measuring time of 1st LS
            startLSTime = omp_get_wtime();

            // Solve the linear system
            tridiagonalSystemSolver(la_y, lb_y, lc_y, c_prime_y, d_prime_y, Ny, LS_b_y, result_y);

            // Update with the result
            for (i = 0; i < Ny; i++)
            {
                index = i * Nx + j;
                RHS[index] = result_y[i];
            }

            // End measuring time of 1st LS
            finishLSTime = omp_get_wtime();
            elapsedTime1stLS += finishLSTime - startLSTime;
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
                index = i * Nx + j;
                actualVm = RHS[index];
                diff_term = diff_coeff * tau * phi_y * (RHS[(lim(i - 1, Ny) * Nx + j)] - 2.0f * actualVm + RHS[(lim(i + 1, Ny) * Nx + j)]);
                LS_b_x[j] = actualVm + diff_term + 0.5f * partRHS[index];
            }

            // Start measuring time of 2nd LS
            startLSTime = omp_get_wtime();

            // Solve the linear system
            tridiagonalSystemSolver(la_x, lb_x, lc_x, c_prime_x, d_prime_x, Nx, LS_b_x, result_x);

            // Update with the result
            for (j = 0; j < Nx; j++)
            {
                index = i * Nx + j;
                Vm[index] = result_x[j];
            }

            // End measuring time of 2nd LS
            finishLSTime = omp_get_wtime();
            elapsedTime2ndLS += finishLSTime - startLSTime;
        }

        // End of the 2nd part of the time step
        finishTime = omp_get_wtime();
        elapsedTime2ndPart += finishTime - startTime;

        // Save frame if needed
        if (saveFrames)
        {
            startSaveFramesTime = omp_get_wtime();

            // If save frames is true and time step is multiple of frame save rate
            if (timeStepCounter % frameSaveRate == 0)
            {
                // Save frame
                snprintf(file_path, MAX_STRING_SIZE, "%s/Vm_%05d.%s", pathToSaveData, timeStepCounter, file_extension);
                save_function(file_path, Vm, Nx, Ny, delta_x, delta_y);
                SUCCESSMSG("Frame at time %.2f ms saved to %s\n", actualTime, file_path);
            }

            finishSaveFramesTime = omp_get_wtime();
            elapsedSaveFramesTime += finishSaveFramesTime - startSaveFramesTime;
        }

        // Measure velocity if needed
        if (measureVelocity)
        {
            startMeasureVelocityTime = omp_get_wtime();

            // Calculate stim velocity
            if (!stim_velocity_measured)
            {
                if (!aux_stim_velocity_flag)
                {
                    point_potential = Vm[begin_point_index];
                    if (point_potential > 10.0f)
                    {
                        begin_point_time = actualTime;
                        aux_stim_velocity_flag = true;
                    }
                }
                else
                {
                    point_potential = Vm[end_point_index];
                    if (point_potential > 10.0f)
                    {
                        end_point_time = actualTime;
                        stim_velocity = (end_point- begin_point) / (end_point_time - begin_point_time); // cm/ms
                        stim_velocity = stim_velocity * 1000.0f;                                         // cm/s
                        stim_velocity_measured = true;
                        INFOMSG("Stim velocity (measured from %.2f to %.2f cm) is %.4g cm/s\n", begin_point, end_point, stim_velocity);
                    }
                }
            }

            finishMeasureVelocityTime = omp_get_wtime();
            elapsedMeasureVelocityTime += finishMeasureVelocityTime - startMeasureVelocityTime;
        }

        // Update time step counter
        timeStepCounter++;
    }

    finishExecutionTime = omp_get_wtime();
    elapsedExecutionTime += finishExecutionTime - startExecutionTime;

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
void run_OSADI_AFHN(const SimulationConfig *config, Measurement *measurement, const real *time_array, real *Vm, real *W, real *partRHS)
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
    const save_function_t save_function = config->save_function;
    const bool measureVelocity = config->measure_velocity;

    // Variables for time measurement
    real startTime = 0.0f;
    real finishTime = 0.0f;
    real startExecutionTime = 0.0f;
    real finishExecutionTime = 0.0f;
    real elapsedExecutionTime = 0.0f;
    real elapsedTime1stPart = 0.0f;
    real elapsedTime2ndPart = 0.0f;
    real startLSTime = 0.0f;
    real finishLSTime = 0.0f;
    real elapsedTime1stLS = 0.0f;
    real elapsedTime2ndLS = 0.0f;
    real startSaveFramesTime = 0.0f;
    real finishSaveFramesTime = 0.0f;
    real elapsedSaveFramesTime = 0.0f;
    real startMeasureVelocityTime = 0.0f;
    real finishMeasureVelocityTime = 0.0f;
    real elapsedMeasureVelocityTime = 0.0f;

    // Auxiliary variables for the loops
    int timeStepCounter = 0;
    real actualTime = 0.0f;
    int i, j, index, si;

    // Save frames variables
    char file_path[MAX_STRING_SIZE];
    char file_extension[4];
    if (strstr(saveFunctionName, "txt") != NULL)
        snprintf(file_extension, 4, "txt");
    else if (strstr(saveFunctionName, "vtk") != NULL)
        snprintf(file_extension, 4, "vtk");

    // Measure velocity variables
    real stim_velocity = 0.0f;
    real point_potential = 0.0f;
    real begin_point = config->Lx / 3.0f;
    real end_point = 2.0f * begin_point;
    int begin_point_index = round(begin_point / delta_x) + 1;
    int end_point_index = round(end_point / delta_x) + 1;
    real begin_point_time = 0.0f;
    real end_point_time = 0.0f;
    bool aux_stim_velocity_flag = false;
    bool stim_velocity_measured = false;

    // Auxiliary variables for the operations
    real stim = 0.0f;
    real diff_term = 0.0f;
    real actualVm;
    real actualW;
    real RHS_Vm_term;

    // Calculate coefficients for the ADI method
    const real phi_x = delta_t / (delta_x * delta_x);
    const real phi_y = delta_t / (delta_y * delta_y);
    const real chi = 1.0e3f;    // cm^-1
    const real Cm = 1.0e-3f;    // mF * cm^-2
    const real diff_coeff = sigma / (Cm * chi);

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
                index = i * Nx + j;

                // Calculate the explicit part of the RHS, including the diffusion term in both directions
                actualVm = Vm[index];
                actualW = W[index];
                RHS_Vm_term = (G * actualVm * (1.0f - (actualVm / vth)) * (1.0f - (actualVm / vp))) + (eta1 * actualVm * actualW);

                // Stimulation
                stim = 0.0f;

#pragma unroll
                for (si = 0; si < numberOfStimuli; si++)
                {
                    if (actualTime >= stimuli[si].begin_time && actualTime <= stimuli[si].begin_time + stimuli[si].duration && j >= stimuli[si].x_discretized.min && j <= stimuli[si].x_discretized.max && i >= stimuli[si].y_discretized.min && i <= stimuli[si].y_discretized.max)
                    {
                        stim = stimuli[si].amplitude;
                        break;
                    }
                }

                // Calculate part of the RHS of the following linear systems with Forward Euler
                partRHS[index] = delta_t * (stim - RHS_Vm_term);

                // Update state variables
                W[index] = actualW + delta_t * (eta2 * ((actualVm / vp) - (eta3 * actualW))); // with Forward Euler -> Wn+1 = Wn + dt*R(Vmn, Wn)
            }
        }

        // End of the 1st part of the time step
        finishTime = omp_get_wtime();
        elapsedTime1stPart += finishTime - startTime;

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
                index = i * Nx + j;
                LS_b_y[i] = Vm[index] + 0.5f * partRHS[index];
            }

            // Start measuring time of 1st LS
            startLSTime = omp_get_wtime();

            // Solve the linear system
            tridiagonalSystemSolver(la_y, lb_y, lc_y, c_prime_y, d_prime_y, Ny, LS_b_y, result_y);

            // Update with the result
            for (i = 0; i < Ny; i++)
            {
                index = i * Nx + j;
                Vm[index] = result_y[i];
            }

            // End measuring time of 1st LS
            finishLSTime = omp_get_wtime();
            elapsedTime1stLS += finishLSTime - startLSTime;
        }

        // ================================================!
        //  Calculate Vm at n+1 -> Result goes to Vm         !
        // ================================================!
        for (i = 0; i < Ny; i++)
        {
            // Calculate the RHS of the linear system
            for (j = 0; j < Nx; j++)
            {
                index = i * Nx + j;
                LS_b_x[j] = Vm[index] + 0.5f * partRHS[index];
            }

            // Start measuring time of 2nd LS
            startLSTime = omp_get_wtime();

            // Solve the linear system
            tridiagonalSystemSolver(la_x, lb_x, lc_x, c_prime_x, d_prime_x, Nx, LS_b_x, result_x);

            // Update with the result
            for (j = 0; j < Nx; j++)
            {
                index = i * Nx + j;
                Vm[index] = result_x[j];
            }

            // End measuring time of 2nd LS
            finishLSTime = omp_get_wtime();
            elapsedTime2ndLS += finishLSTime - startLSTime;
        }

        // End of the 2nd part of the time step
        finishTime = omp_get_wtime();
        elapsedTime2ndPart += finishTime - startTime;

        // Save frame if needed
        if (saveFrames)
        {
            startSaveFramesTime = omp_get_wtime();

            // If save frames is true and time step is multiple of frame save rate
            if (timeStepCounter % frameSaveRate == 0)
            {
                // Save frame
                snprintf(file_path, MAX_STRING_SIZE, "%s/Vm_%05d.%s", pathToSaveData, timeStepCounter, file_extension);
                save_function(file_path, Vm, Nx, Ny, delta_x, delta_y);
                SUCCESSMSG("Frame at time %.2f ms saved to %s\n", actualTime, file_path);
            }

            finishSaveFramesTime = omp_get_wtime();
            elapsedSaveFramesTime += finishSaveFramesTime - startSaveFramesTime;
        }

        // Measure velocity if needed
        if (measureVelocity)
        {
            startMeasureVelocityTime = omp_get_wtime();

            // Calculate stim velocity
            if (!stim_velocity_measured)
            {
                if (!aux_stim_velocity_flag)
                {
                    point_potential = Vm[begin_point_index];
                    if (point_potential > 10.0f)
                    {
                        begin_point_time = actualTime;
                        aux_stim_velocity_flag = true;
                    }
                }
                else
                {
                    point_potential = Vm[end_point_index];
                    if (point_potential > 10.0f)
                    {
                        end_point_time = actualTime;
                        stim_velocity = (end_point- begin_point) / (end_point_time - begin_point_time); // cm/ms
                        stim_velocity = stim_velocity * 1000.0f;                                         // cm/s
                        stim_velocity_measured = true;
                        INFOMSG("Stim velocity (measured from %.2f to %.2f cm) is %.4g cm/s\n", begin_point, end_point, stim_velocity);
                    }
                }
            }

            finishMeasureVelocityTime = omp_get_wtime();
            elapsedMeasureVelocityTime += finishMeasureVelocityTime - startMeasureVelocityTime;
        }

        // Update time step counter
        timeStepCounter++;
    }

    finishExecutionTime = omp_get_wtime();
    elapsedExecutionTime += finishExecutionTime - startExecutionTime;

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
void run_FE_AFHN(const SimulationConfig *config, Measurement *measurement, const real *time_array, real *Vm, real *W, real *partRHS)
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
    const save_function_t save_function = config->save_function;
    const bool measureVelocity = config->measure_velocity;

    // Variables for time measurement
    real startTime = 0.0f;
    real finishTime = 0.0f;
    real startExecutionTime = 0.0f;
    real finishExecutionTime = 0.0f;
    real elapsedExecutionTime = 0.0f;
    real elapsedTime1stPart = 0.0f;
    real elapsedTime2ndPart = 0.0f;
    real startSaveFramesTime = 0.0f;
    real finishSaveFramesTime = 0.0f;
    real elapsedSaveFramesTime = 0.0f;
    real startMeasureVelocityTime = 0.0f;
    real finishMeasureVelocityTime = 0.0f;
    real elapsedMeasureVelocityTime = 0.0f;

    // Auxiliary variables for the loops
    int timeStepCounter = 0;
    real actualTime = 0.0f;
    int i, j, index, si;

    // Save frames variables
    char file_path[MAX_STRING_SIZE];
    char file_extension[4];
    if (strstr(saveFunctionName, "txt") != NULL)
        snprintf(file_extension, 4, "txt");
    else if (strstr(saveFunctionName, "vtk") != NULL)
        snprintf(file_extension, 4, "vtk");

    // Measure velocity variables
    real stim_velocity = 0.0f;
    real point_potential = 0.0f;
    real begin_point = config->Lx / 3.0f;
    real end_point = 2.0f * begin_point;
    int begin_point_index = round(begin_point / delta_x) + 1;
    int end_point_index = round(end_point / delta_x) + 1;
    real begin_point_time = 0.0f;
    real end_point_time = 0.0f;
    bool aux_stim_velocity_flag = false;
    bool stim_velocity_measured = false;

    // Auxiliary variables for the operations
    real stim = 0.0f;
    real diff_term = 0.0f;
    real actualVm;
    real actualW;
    real RHS_Vm_term;

    // Calculate coefficients for the ADI method
    const real phi_x = delta_t / (delta_x * delta_x);
    const real phi_y = delta_t / (delta_y * delta_y);
    const real chi = 1.0e3f;    // cm^-1
    const real Cm = 1.0e-3f;    // mF * cm^-2
    const real diff_coeff = sigma / (Cm * chi);

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
                index = i * Nx + j;

                // Calculate the explicit part of the RHS, including the diffusion term in both directions
                actualVm = Vm[index];
                actualW = W[index];
                RHS_Vm_term = (G * actualVm * (1.0f - (actualVm / vth)) * (1.0f - (actualVm / vp))) + (eta1 * actualVm * actualW);

                // Stimulation
                stim = 0.0f;

#pragma unroll
                for (si = 0; si < numberOfStimuli; si++)
                {
                    if (actualTime >= stimuli[si].begin_time && actualTime <= stimuli[si].begin_time + stimuli[si].duration && j >= stimuli[si].x_discretized.min && j <= stimuli[si].x_discretized.max && i >= stimuli[si].y_discretized.min && i <= stimuli[si].y_discretized.max)
                    {
                        stim = stimuli[si].amplitude;
                        break;
                    }
                }

                // Update variables explicitly
                diff_term = diff_coeff * (phi_x * (Vm[i * Nx + (lim(j - 1, Nx))] - 2.0f * actualVm + Vm[i * Nx + lim(j + 1, Nx)]) + phi_y * (Vm[(lim(i - 1, Ny) * Nx + j)] - 2.0f * actualVm + Vm[(lim(i + 1, Ny) * Nx + j)]));
                partRHS[index] = actualVm + diff_term + delta_t * (stim - RHS_Vm_term);

                W[index] = actualW + delta_t * (eta2 * ((actualVm / vp) - (eta3 * actualW)));
            }
        }

        // End of the 1st part of the time step
        finishTime = omp_get_wtime();
        elapsedTime1stPart += finishTime - startTime;

        // Start measuring time of 2nd part
        startTime = omp_get_wtime();

        // ==================!
        //  Update Vm         !
        // ==================!
        for (i = 0; i < Ny; i++)
            for (j = 0; j < Nx; j++)
            {
                index = i * Nx + j;
                Vm[index] = partRHS[index];
            }

        // End of the 2nd part of the time step
        finishTime = omp_get_wtime();
        elapsedTime2ndPart += finishTime - startTime;

        // Save frame if needed
        if (saveFrames)
        {
            startSaveFramesTime = omp_get_wtime();

            // If save frames is true and time step is multiple of frame save rate
            if (timeStepCounter % frameSaveRate == 0)
            {
                // Save frame
                snprintf(file_path, MAX_STRING_SIZE, "%s/Vm_%05d.%s", pathToSaveData, timeStepCounter, file_extension);
                save_function(file_path, Vm, Nx, Ny, delta_x, delta_y);
                SUCCESSMSG("Frame at time %.2f ms saved to %s\n", actualTime, file_path);
            }

            finishSaveFramesTime = omp_get_wtime();
            elapsedSaveFramesTime += finishSaveFramesTime - startSaveFramesTime;
        }

        // Measure velocity if needed
        if (measureVelocity)
        {
            startMeasureVelocityTime = omp_get_wtime();

            // Calculate stim velocity
            if (!stim_velocity_measured)
            {
                if (!aux_stim_velocity_flag)
                {
                    point_potential = Vm[begin_point_index];
                    if (point_potential > 10.0f)
                    {
                        begin_point_time = actualTime;
                        aux_stim_velocity_flag = true;
                    }
                }
                else
                {
                    point_potential = Vm[end_point_index];
                    if (point_potential > 10.0f)
                    {
                        end_point_time = actualTime;
                        stim_velocity = (end_point- begin_point) / (end_point_time - begin_point_time); // cm/ms
                        stim_velocity = stim_velocity * 1000.0f;                                         // cm/s
                        stim_velocity_measured = true;
                        INFOMSG("Stim velocity (measured from %.2f to %.2f cm) is %.4g cm/s\n", begin_point, end_point, stim_velocity);
                    }
                }
            }

            finishMeasureVelocityTime = omp_get_wtime();
            elapsedMeasureVelocityTime += finishMeasureVelocityTime - startMeasureVelocityTime;
        }

        // Update time step counter
        timeStepCounter++;
    }

    finishExecutionTime = omp_get_wtime();
    elapsedExecutionTime += finishExecutionTime - startExecutionTime;

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
