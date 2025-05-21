#include "../numerical_methods.h"
#include "../numerical_methods_helpers.h"
#include "../../cell_models/cell_models.h"

static __global__ void computeReactionApproxAndUpdateSV(const int Nx, const int Ny, const real delta_t, const real phi_x, const real phi_y,
                                                               const real diff_coeff, const real actualTime, const int num_active_stimuli,
                                                               const Stimulus *d_active_stimuli, const real *d_Vm, real *d_sV, real *d_partRHS, const CellModel cell_model)
{
    // Obtain the thread index
    const int i = blockIdx.y * blockDim.y + threadIdx.y;
    const int j = blockIdx.x * blockDim.x + threadIdx.x;
    const int idx = i * Nx + j;

    if (i < Ny && j < Nx)
    {
        // Declare auxiliary arrays
        real d_actualsV[MAX_NSV];
        real d_sVtilde[MAX_NSV];

        // Calculate the explicit part of the RHS, including the diffusion term in both directions
        real actualVm = d_Vm[idx];
        select_get_actual_sV(cell_model, d_actualsV, d_sV, idx);
        
        // Stimulation
        real stim = get_stimulus_value(actualTime, i, j, d_active_stimuli, num_active_stimuli);

        // Calculate aproximation with RK2 -> Vmn+1/2 = Vmn + 0.5*diffusion + 0.5*dt*R(Vmn, Wn)
        real diff_term = compute_diffusion_term(d_Vm, i, j, Nx, Ny, diff_coeff, phi_x, phi_y);
        real Vmtilde = actualVm + 0.5f * diff_term + (0.5f * delta_t * (stim - select_compute_dVmdt(cell_model, actualVm, d_actualsV)));

        // Calculate approximation for state variables and prepare part of the RHS of the following linear systems
        select_update_sVtilde(cell_model, d_sVtilde, actualVm, d_actualsV, 0.5f * delta_t);
        d_partRHS[idx] = delta_t * (stim - select_compute_dVmdt(cell_model, Vmtilde, d_sVtilde));

        // Update state variables
        select_update_sV(cell_model, d_sV, d_actualsV, Vmtilde, d_sVtilde, delta_t, idx);
    }
}

static __global__ void prepareRHS_x(const int Nx, const int Ny, const real coeff,
                                    const real *d_Vm, const real *d_partRHS, real *d_RHS)
{
    // Obtain the thread index
    const int i = blockIdx.y * blockDim.y + threadIdx.y;
    const int j = blockIdx.x * blockDim.x + threadIdx.x;
    const int idx = i * Nx + j;
    const int idx_left = i * Nx + lim(j - 1, Nx);
    const int idx_right = i * Nx + lim(j + 1, Nx);

    if (i < Ny && j < Nx)
    {
        // Calculate the RHS of the linear system with the explicit diffusion term along x
        const real Vm_center = d_Vm[idx];
        const real Vm_left = d_Vm[idx_left];
        const real Vm_right = d_Vm[idx_right];
        const real diff_term = coeff * (Vm_left - 2.0f * Vm_center + Vm_right);
        d_RHS[idx] = Vm_center + diff_term + 0.5f * d_partRHS[idx];
    }
}

static __global__ void prepareRHS_y(const int Nx, const int Ny, const real coeff,
                                    const real *d_Vm, const real *d_partRHS, real *d_RHS)
{
    // Obtain the thread index
    const int i = blockIdx.y * blockDim.y + threadIdx.y;
    const int j = blockIdx.x * blockDim.x + threadIdx.x;
    const int idx = i * Nx + j;
    const int idx_top = lim(i + 1, Ny) * Nx + j;
    const int idx_bottom = lim(i - 1, Ny) * Nx + j;

    if (i < Ny && j < Nx)
    {
        // Calculate the RHS of the linear system with the explicit diffusion term along y
        const real Vm_center = d_Vm[idx];
        const real Vm_top = d_Vm[idx_top];
        const real Vm_bottom = d_Vm[idx_bottom];
        const real diff_term = coeff * (Vm_bottom - 2.0f * Vm_center + Vm_top);
        d_RHS[idx] = Vm_center + diff_term + 0.5f * d_partRHS[idx];
    }
}

static __global__ void parallelThomas_x(const int numSys, const int sysSize, real *d_rhs,
                                        const real *d_la, const real *d_lb, const real *d_lc)
{
    // Obtain the index of the thread - each thread will handle a system
    const int sysIdx = blockIdx.x * blockDim.x + threadIdx.x;
    const int offset = sysIdx * sysSize;

    if (sysIdx < numSys)
    {
        // Local variables
        real c_prime[MAX_SYS_SIZE];
        real d_prime[MAX_SYS_SIZE];
        c_prime[0] = d_lc[0] / d_lb[0];
        d_prime[0] = d_rhs[offset] / d_lb[0];

        for (int i = 1; i < sysSize; i++)
        {
            const real la = d_la[i];
            const real lb = d_lb[i];
            const real lc = d_lc[i];
            const real denom = 1.0f / (lb - c_prime[i - 1] * la);
            if (i < sysSize - 1)
                c_prime[i] = lc * denom;
            d_prime[i] = (d_rhs[offset + i] - d_prime[i - 1] * la) * denom;
        }

        d_rhs[offset + sysSize - 1] = d_prime[sysSize - 1];

        for (int i = sysSize - 2; i >= 0; i--)
            d_rhs[offset + i] = d_prime[i] - c_prime[i] * d_rhs[offset + i + 1];
    }
}

static __global__ void parallelThomas_y(const int numSys, const int sysSize, real *d_rhs,
                                        const real *d_la, const real *d_lb, const real *d_lc)
{
    // Obtain the index of the thread - each thread will handle a system
    const int sysIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (sysIdx < numSys)
    {
        // Local variables
        real c_prime[MAX_SYS_SIZE];
        real d_prime[MAX_SYS_SIZE];
        c_prime[0] = d_lc[0] / d_lb[0];
        d_prime[0] = d_rhs[sysIdx] / d_lb[0];

        for (int i = 1; i < sysSize; i++)
        {
            const real la = d_la[i];
            const real lb = d_lb[i];
            const real lc = d_lc[i];
            const real denom = 1.0f / (lb - c_prime[i - 1] * la);
            if (i < sysSize - 1)
                c_prime[i] = lc * denom;
            d_prime[i] = (d_rhs[sysIdx + numSys * i] - d_prime[i - 1] * la) * denom;
        }

        d_rhs[sysIdx + numSys * (sysSize - 1)] = d_prime[sysSize - 1];

        for (int i = sysSize - 2; i >= 0; i--)
            d_rhs[sysIdx + numSys * i] = d_prime[i] - c_prime[i] * d_rhs[sysIdx + numSys * (i + 1)];
    }
}

// TODO: validate this function
// static __global__ void parallelThomas_large(const int numSys, const int sysSize, real *d_rhs, const real *d_la,
//                                             const real *d_lb, const real *d_lc)
// {
//     const int sysIdx = blockIdx.x * blockDim.x + threadIdx.x;
//     if (sysIdx >= numSys)
//         return;

//     const int offset = sysIdx * sysSize;

//     // Forward pass (optimized for large systems)
//     const real c_prime_0 = d_lc[0] / d_lb[0];
//     real c_prime_prev = c_prime_0;
//     real d_prime_prev = d_rhs[offset] / d_lb[0];
//     d_rhs[offset] = d_prime_prev; // Reuse d_rhs for temporary storage

//     for (int i = 1; i < sysSize; ++i)
//     {
//         const real la = d_la[i];
//         const real lb = d_lb[i];
//         const real lc = d_lc[i];
//         const real rhs = d_rhs[offset + i];
//         const real denom = 1.0f / (lb - c_prime_prev * la);
//         const real new_c_prime = (i < sysSize - 1) ? lc * denom : 0.0f;
//         const real new_d_prime = (rhs - d_prime_prev * la) * denom;

//         d_rhs[offset + i] = new_d_prime;
//         c_prime_prev = new_c_prime;
//         d_prime_prev = new_d_prime;
//     }

//     // Backward pass
//     for (int i = sysSize - 2; i >= 0; --i)
//     {
//         d_rhs[offset + i] -= d_rhs[offset + i + 1] * c_prime_prev;
//         // Update c_prime_prev for next iteration
//         if (i > 0)
//         {
//             c_prime_prev = d_lc[i - 1] / (d_lb[i - 1] - d_la[i - 1] * c_prime_prev);
//         }
//         else
//         {
//             c_prime_prev = c_prime_0;
//         }
//     }
// }

void runSSIADI_CUDA(const SimulationConfig *config, Measurement *measurement, const real *time_array,
                    const CellModelSolver *cell_model_solver, real *Vm, real *sV)
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

    const CellModel cell_model = config->cell_model;
    const bool saveFrames = config->save_frames;
    const int frameSaveRate = config->frame_save_rate;
    const char *pathToSaveData = config->output_dir;
    const char *file_extension = config->file_extension;
    const save_function_t save_function = config->save_function;
    const bool measureVelocity = config->measure_velocity;
    
    // Get the solver functions
    const real activation_threshold = cell_model_solver->activation_thershold;
    const compute_diffusion_coefficient_t compute_diffusion_coefficient = cell_model_solver->compute_diffusion_coefficient;
    
    // Measure velocity variables
    real stim_velocity, t0, t1;
    const real x0 = config->Lx / 3.0f;
    const real x1 = 2.0f * x0;
    const int idx_x0 = round(x0 / delta_x) + 1;
    const int idx_x1 = round(x1 / delta_x) + 1;
    bool aux_stim_velocity_flag = false;
    bool stim_velocity_measured = false;

    // Auxiliary variables for the loops
    int timeStepCounter = 0;
    real actualTime = 0.0f;
    
    // Create device variables, allocate memory on device, and copy data
    real *d_Vm, *d_sV;
    Stimulus *d_stimuli;
    const int total_points = Nx * Ny;

    CUDA_CALL(cudaMalloc(&d_Vm, total_points * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_sV, total_points * cell_model_solver->n_state_vars * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_stimuli, numberOfStimuli * sizeof(Stimulus)));

    CUDA_CALL(cudaMemcpy(d_Vm, Vm, total_points * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_sV, sV, total_points * cell_model_solver->n_state_vars * sizeof(real), cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(d_stimuli, stimuli, numberOfStimuli * sizeof(Stimulus), cudaMemcpyHostToDevice));
    
    // Auxiliary variables for the operations
    real *d_partRHS, *d_RHS;
    CUDA_CALL(cudaMalloc(&d_partRHS, total_points * sizeof(real)));
    CUDA_CALL(cudaMalloc(&d_RHS, total_points * sizeof(real)));
    
    // Calculate coefficients for the ADI method
    const real phi_x = delta_t / (delta_x * delta_x);
    const real phi_y = delta_t / (delta_y * delta_y);
    const real diff_coeff = compute_diffusion_coefficient(sigma);
    const real tau = 0.5f; // Used for the explicit RHS of ADI
    
    // Auxiliary arrays for Thomas algorithm
    real *la_x = (real *)malloc(Nx * sizeof(real)); // subdiagonal
    real *lb_x = (real *)malloc(Nx * sizeof(real)); // diagonal
    real *lc_x = (real *)malloc(Nx * sizeof(real)); // superdiagonal
    real *la_y = (real *)malloc(Ny * sizeof(real)); // subdiagonal
    real *lb_y = (real *)malloc(Ny * sizeof(real)); // diagonal
    real *lc_y = (real *)malloc(Ny * sizeof(real)); // superdiagonal
    
    populateDiagonalThomasAlgorithm(la_x, lb_x, lc_x, Nx, 0.5f * phi_x * diff_coeff);
    populateDiagonalThomasAlgorithm(la_y, lb_y, lc_y, Ny, 0.5f * phi_y * diff_coeff);
    
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
    
    // CUDA grid and block allocation
    // Device properties
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);

    // Number of SMs and minimum number of blocks to maximize the parallelism
    const int numSMs = prop.multiProcessorCount;
    const int minBlocks = 2 * numSMs;

    // Print information
    SIMPLEMSG("");
    INFOMSG("Device name: %s (%d SMs)\n", prop.name, numSMs);

    // Calculate the number of blocks and threads for the full domain kernels
    dim3 fullDomainBlockSize(FULL_DOMAIN_BLOCK_SIZE_X, FULL_DOMAIN_BLOCK_SIZE_Y);
    dim3 fullDomainGridSize((Nx + fullDomainBlockSize.x - 1) / fullDomainBlockSize.x, (Ny + fullDomainBlockSize.y - 1) / fullDomainBlockSize.y);

    // Adjust the number of blocks
    if (fullDomainGridSize.x * fullDomainGridSize.y < minBlocks)
        fullDomainGridSize.x = (minBlocks + fullDomainGridSize.y - 1) / fullDomainGridSize.y;
    
    // Print information
    SIMPLEMSG("");
    INFOMSG("For full domain kernels:\n");
    INFOMSG("Block size: %d x %d threads (total %d threads per block)\n", fullDomainBlockSize.x, fullDomainBlockSize.y, fullDomainBlockSize.x * fullDomainBlockSize.y);
    INFOMSG("Grid size: %d x %d blocks (total %d blocks, total %d threads)\n", fullDomainGridSize.x, fullDomainGridSize.y, fullDomainGridSize.x * fullDomainGridSize.y, fullDomainGridSize.x * fullDomainGridSize.y * fullDomainBlockSize.x * fullDomainBlockSize.y);

    // Calculate the number of blocks and threads for individual directions of ADI that will be used in Thomas kernel
    const int gridSize_x = (Nx + THOMAS_KERNEL_BLOCK_SIZE - 1) / THOMAS_KERNEL_BLOCK_SIZE;
    const int gridSize_y = (Ny + THOMAS_KERNEL_BLOCK_SIZE - 1) / THOMAS_KERNEL_BLOCK_SIZE;

    // Print information
    SIMPLEMSG("");
    INFOMSG("For Thomas kernel:\n");
    INFOMSG("Grid size for x: %d blocks (%d threads per block, total %d threads)\n", gridSize_x, THOMAS_KERNEL_BLOCK_SIZE, gridSize_x * THOMAS_KERNEL_BLOCK_SIZE);
    INFOMSG("Grid size for y: %d blocks (%d threads per block, total %d threads)\n", gridSize_y, THOMAS_KERNEL_BLOCK_SIZE, gridSize_y * THOMAS_KERNEL_BLOCK_SIZE);

    // Variables for time measurement
    real startTime = 0.0f;
    real startExecutionTime = 0.0f;
    real elapsedExecutionTime = 0.0f;
    real elapsedTime1stPart = 0.0f;
    real elapsedTime2ndPart = 0.0f;
    real elapsedSaveFramesTime = 0.0f;
    real elapsedMeasureVelocityTime = 0.0f;
    real startLSTime = 0.0f;
    real elapsedTime1stLS = 0.0f;
    real elapsedTime2ndLS = 0.0f;

    SIMPLEMSG("");
    INFOMSG("Starting simulation with SSIADI (CUDA)...\n");

    // Main time loop
    startExecutionTime = omp_get_wtime();

    while (timeStepCounter < M)
    {
        // Get time step
        actualTime = time_array[timeStepCounter];

        // ================================================!
        //  Calculate Approxs. and Update ODEs             !
        // ================================================!
        startTime = omp_get_wtime();

        // Launch kernel to compute the reaction term and update state variables
        computeReactionApproxAndUpdateSV<<<fullDomainGridSize, fullDomainBlockSize>>>(Nx, Ny, delta_t, phi_x, phi_y, diff_coeff, actualTime,
                                                                                 numberOfStimuli, d_stimuli, d_Vm, d_sV, d_partRHS, cell_model);
        CUDA_CALL(cudaDeviceSynchronize());

        elapsedTime1stPart += omp_get_wtime() - startTime;

        // ================================================!
        //  Calculate Vm at n+1/2 -> Result goes to RHS    !
        //  diffusion implicit in y and explicit in x      !
        // ================================================!
        startTime = omp_get_wtime();

        prepareRHS_x<<<fullDomainGridSize, fullDomainBlockSize>>>(Nx, Ny, phi_x * diff_coeff * tau, d_Vm, d_partRHS, d_RHS);
        CUDA_CALL(cudaDeviceSynchronize());

        startLSTime = omp_get_wtime();

        parallelThomas_y<<<gridSize_x, THOMAS_KERNEL_BLOCK_SIZE>>>(Nx, Ny, d_RHS, d_la_y, d_lb_y, d_lc_y);
        CUDA_CALL(cudaDeviceSynchronize());

        elapsedTime1stLS += omp_get_wtime() - startLSTime;

        // ================================================!
        //  Calculate Vm at n+1 -> Result goes to Vm       !
        //  diffusion implicit in x and explicit in y      !
        // ================================================!
        prepareRHS_y<<<fullDomainGridSize, fullDomainBlockSize>>>(Nx, Ny, phi_y * diff_coeff * tau, d_RHS, d_partRHS, d_Vm);
        CUDA_CALL(cudaDeviceSynchronize());

        startLSTime = omp_get_wtime();

        parallelThomas_x<<<gridSize_y, THOMAS_KERNEL_BLOCK_SIZE>>>(Ny, Nx, d_Vm, d_la_x, d_lb_x, d_lc_x);
        CUDA_CALL(cudaDeviceSynchronize());

        elapsedTime2ndLS += omp_get_wtime() - startLSTime;

        elapsedTime2ndPart += omp_get_wtime() - startTime;

        // Save frame if needed
        if (saveFrames && (timeStepCounter % frameSaveRate == 0))
        {
            startTime = omp_get_wtime();
            CUDA_CALL(cudaMemcpy(Vm, d_Vm, total_points * sizeof(real), cudaMemcpyDeviceToHost));
            handle_frame_saving(pathToSaveData, file_extension, save_function, timeStepCounter, Vm, Nx, Ny, delta_x, delta_y, actualTime);
            elapsedSaveFramesTime += omp_get_wtime() - startTime;
        }

        // Measure velocity if needed
        if (measureVelocity && !stim_velocity_measured)
        {
            startTime = omp_get_wtime();

            // Copy only the Vm values needed for velocity measurement
            real Vmidx_x0, Vmidx_x1;
            CUDA_CALL(cudaMemcpy(&Vmidx_x0, &d_Vm[idx_x0], sizeof(real), cudaMemcpyDeviceToHost));
            CUDA_CALL(cudaMemcpy(&Vmidx_x1, &d_Vm[idx_x1], sizeof(real), cudaMemcpyDeviceToHost));
            handle_velocity_measurement(Vmidx_x0, Vmidx_x1, &t0, &t1, activation_threshold, &aux_stim_velocity_flag, &stim_velocity_measured, actualTime, x0, x1, &stim_velocity);
            
            elapsedMeasureVelocityTime += omp_get_wtime() - startTime;
        }
        
        // Update time step counter
        timeStepCounter++;
    }

    elapsedExecutionTime = omp_get_wtime() - startExecutionTime;

    INFOMSG("Simulation done!\n");
    SIMPLEMSG("");

    // Update measurement structure
    measurement->elapsedExecutionTime = elapsedExecutionTime;
    measurement->elapsedTime1stPart = elapsedTime1stPart;
    measurement->elapsedTime2ndPart = elapsedTime2ndPart;
    measurement->elapsedTime1stLS = elapsedTime1stLS;
    measurement->elapsedTime2ndLS = elapsedTime2ndLS;
    measurement->elapsedSaveFramesTime = elapsedSaveFramesTime;
    measurement->elapsedMeasureVelocityTime = elapsedMeasureVelocityTime;
    measurement->stimVelocity = stim_velocity;

    // Copy results back to host
    CUDA_CALL(cudaMemcpy(Vm, d_Vm, total_points * sizeof(real), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(sV, d_sV, total_points * cell_model_solver->n_state_vars * sizeof(real), cudaMemcpyDeviceToHost));

    // Free allocated memory
    CUDA_CALL(cudaFree(d_Vm));
    CUDA_CALL(cudaFree(d_sV));
    CUDA_CALL(cudaFree(d_stimuli));
    CUDA_CALL(cudaFree(d_partRHS));
    CUDA_CALL(cudaFree(d_RHS));
    CUDA_CALL(cudaFree(d_la_x));
    CUDA_CALL(cudaFree(d_lb_x));
    CUDA_CALL(cudaFree(d_lc_x));
    CUDA_CALL(cudaFree(d_la_y));
    CUDA_CALL(cudaFree(d_lb_y));
    CUDA_CALL(cudaFree(d_lc_y));
}
