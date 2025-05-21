#include "tt2.h"

#ifndef USE_CUDA

void solveMonodomainTT2(const SimulationConfig *config, Measurement *measurement, const real *time_array)
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
    // initializeWithInitialConditionAFHN(Nx, Ny, Vm, W);

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
    // if (config->method == METHOD_SSIADI)
    // {
    //     // Run the simulation using the SSIADI method
    //     run_SSIADI_AFHN(config, measurement, time_array, Vm, W, partRHS);
    // }
    // else if (config->method == METHOD_OSADI)
    // {
    //     // Run the simulation using the OSADI method
    //     run_OSADI_AFHN(config, measurement, time_array, Vm, W, partRHS);
    // }
    // else if (config->method == METHOD_FE)
    // {
    //     // Run the simulation using the Forward Euler method
    //     run_FE_AFHN(config, measurement, time_array, Vm, W, partRHS);
    // }
    // else
    // {
    //     fprintf(stderr, "Error: Unsupported method for AFHN model.\n");
    //     return;
    // }

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
        snprintf(file_path, MAX_STRING_SIZE, "%s/frames/Vm_%05d.%s", config->output_dir, M, file_extension);
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

#endif // USE_CUDA