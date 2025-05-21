#include "../numerical_methods.h"
#include "../numerical_methods_helpers.h"
#include "../../cell_models/cell_models.h"

static __global__ void computeReactionAndUpdateSV(const int Nx, const int Ny, const real delta_t, const real actualTime, const int num_active_stimuli,
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

        // Calculate the explicit part of the RHS, including the diffusion term in both directions
        real actualVm = d_Vm[idx];
        select_get_actual_sV(cell_model, d_actualsV, d_sV, idx);
        
        // Stimulation
        real stim = get_stimulus_value(actualTime, i, j, d_active_stimuli, num_active_stimuli);

        // Calculate part of the RHS of the following linear systems with Forward Euler
        d_partRHS[idx] = delta_t * (stim - select_compute_dVmdt(cell_model, actualVm, d_actualsV));

        // Update state variables
        select_update_sV(cell_model, d_sV, d_actualsV, actualVm, d_actualsV, delta_t, idx);
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

static __global__ void prepareRHS(const int Nx, const int Ny, real *d_Vm, const real *d_partRHS)
{
    // Obtain the thread index
    const int i = blockIdx.y * blockDim.y + threadIdx.y;
    const int j = blockIdx.x * blockDim.x + threadIdx.x;
    const int idx = i * Nx + j;

    if (i < Ny && j < Nx)
        d_Vm[idx] = d_Vm[idx] + 0.5f * d_partRHS[idx];
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

void runOSADI_CUDA(const SimulationConfig *config, Measurement *measurement, const real *time_array,
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
    real *d_partRHS;
    CUDA_CALL(cudaMalloc(&d_partRHS, total_points * sizeof(real)));
    
    // Calculate coefficients for the ADI method
    const real phi_x = delta_t / (delta_x * delta_x);
    const real phi_y = delta_t / (delta_y * delta_y);
    const real diff_coeff = compute_diffusion_coefficient(sigma);
    
    // Auxiliary arrays for Thomas algorithm
    real *la_x = (real *)malloc(Nx * sizeof(real)); // subdiagonal
    real *lb_x = (real *)malloc(Nx * sizeof(real)); // diagonal
    real *lc_x = (real *)malloc(Nx * sizeof(real)); // superdiagonal
    real *la_y = (real *)malloc(Ny * sizeof(real)); // subdiagonal
    real *lb_y = (real *)malloc(Ny * sizeof(real)); // diagonal
    real *lc_y = (real *)malloc(Ny * sizeof(real)); // superdiagonal
    
    populateDiagonalThomasAlgorithm(la_x, lb_x, lc_x, Nx, phi_x * diff_coeff);
    populateDiagonalThomasAlgorithm(la_y, lb_y, lc_y, Ny, phi_y * diff_coeff);
    
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
    INFOMSG("Starting simulation with OSADI (CUDA)...\n");

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
        computeReactionAndUpdateSV<<<fullDomainGridSize, fullDomainBlockSize>>>(Nx, Ny, delta_t, actualTime, numberOfStimuli, d_stimuli, d_Vm, d_sV, d_partRHS, cell_model);
        CUDA_CALL(cudaDeviceSynchronize());

        elapsedTime1stPart += omp_get_wtime() - startTime;

        // ================================================!
        //  Calculate Vm at n+1/2 -> Result goes to Vm     !
        // ================================================!
        startTime = omp_get_wtime();

        prepareRHS<<<fullDomainGridSize, fullDomainBlockSize>>>(Nx, Ny, d_Vm, d_partRHS);
        CUDA_CALL(cudaDeviceSynchronize());

        startLSTime = omp_get_wtime();

        parallelThomas_y<<<gridSize_x, THOMAS_KERNEL_BLOCK_SIZE>>>(Nx, Ny, d_Vm, d_la_y, d_lb_y, d_lc_y);
        CUDA_CALL(cudaDeviceSynchronize());

        elapsedTime1stLS += omp_get_wtime() - startLSTime;

        // ================================================!
        //  Calculate Vm at n+1 -> Result goes to Vm       !
        // ================================================!
        prepareRHS<<<fullDomainGridSize, fullDomainBlockSize>>>(Nx, Ny, d_Vm, d_partRHS);
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
    CUDA_CALL(cudaFree(d_la_x));
    CUDA_CALL(cudaFree(d_lb_x));
    CUDA_CALL(cudaFree(d_lc_x));
    CUDA_CALL(cudaFree(d_la_y));
    CUDA_CALL(cudaFree(d_lb_y));
    CUDA_CALL(cudaFree(d_lc_y));
}
