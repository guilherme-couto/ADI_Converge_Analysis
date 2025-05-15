#include "mv.h"

#ifndef USE_CUDA

// Parameters - Based on Minimal Ventricular model
// Model definition https://www.sciencedirect.com/science/article/pii/S0022519308001690?via%3Dihub

// const real Dtilde = 0.65f * 1e-3; // cm^2/s

#ifdef EPI

const real u_o = 0.0f;
const real u_u = 1.55f;
const real theta_v = 0.3f;
const real theta_w = 0.13f;
const real theta_vminus = 0.006f;
const real theta_o = 0.006f;
const real tau_v1minus = 60.0f;
const real tau_v2minus = 1150.0f;
const real tau_vplus = 1.4506f;
const real tau_w1minus = 60.0f;
const real tau_w2minus = 15.0f;
const real k_wminus = 65.0f;
const real u_wminus = 0.03f;
const real tau_wplus = 200.0f;
const real tau_fi = 0.11f;
const real tau_o1 = 400.0f;
const real tau_o2 = 6.0f;
const real tau_so1 = 30.0181f;
const real tau_so2 = 0.9957f;
const real k_so = 2.0458f;
const real u_so = 0.65f;
const real tau_s1 = 2.7342f;
const real tau_s2 = 16.0f;
const real k_s = 2.0994f;
const real u_s = 0.9087f;
const real tau_si = 1.8875f;
const real tau_winf = 0.07f;
const real w_infstar = 0.94f;

#endif // EPI

#ifdef ENDO

const real u_o = 0.0f;
const real u_u = 1.56f;
const real theta_v = 0.3f;
const real theta_w = 0.13f;
const real theta_vminus = 0.2f;
const real theta_o = 0.006f;
const real tau_v1minus = 75.0f;
const real tau_v2minus = 10.0f;
const real tau_vplus = 1.4506f;
const real tau_w1minus = 6.0f;
const real tau_w2minus = 140.0f;
const real k_wminus = 200.0f;
const real u_wminus = 0.016f;
const real tau_wplus = 280.0f;
const real tau_fi = 0.1f;
const real tau_o1 = 470.0f;
const real tau_o2 = 6.0f;
const real tau_so1 = 40.0f;
const real tau_so2 = 1.2f;
const real k_so = 2.0f;
const real u_so = 0.65f;
const real tau_s1 = 2.7342f;
const real tau_s2 = 2.0f;
const real k_s = 2.0994f;
const real u_s = 0.9087f;
const real tau_si = 2.9013f;
const real tau_winf = 0.0273f;
const real w_infstar = 0.78f;

#endif // ENDO

#ifdef MCELL

const real u_o = 0.0f;
const real u_u = 1.61f;
const real theta_v = 0.3f;
const real theta_w = 0.13f;
const real theta_vminus = 0.1f;
const real theta_o = 0.005f;
const real tau_v1minus = 80.0f;
const real tau_v2minus = 1.4506f;
const real tau_vplus = 1.4506f;
const real tau_w1minus = 70.0f;
const real tau_w2minus = 8.0f;
const real k_wminus = 200.0f;
const real u_wminus = 0.016f;
const real tau_wplus = 280.0f;
const real tau_fi = 0.078f;
const real tau_o1 = 410.0f;
const real tau_o2 = 7.0f;
const real tau_so1 = 91.0f;
const real tau_so2 = 0.8f;
const real k_so = 2.1f;
const real u_so = 0.6f;
const real tau_s1 = 2.7342f;
const real tau_s2 = 4.0f;
const real k_s = 2.0994f;
const real u_s = 0.9087f;
const real tau_si = 3.3849f;
const real tau_winf = 0.01f;
const real w_infstar = 0.5f;

#endif // MCELL

#ifdef PB

const real u_o = 0.0f;
const real u_u = 1.45f;
const real theta_v = 0.35f;
const real theta_w = 0.13f;
const real theta_vminus = 0.175f;
const real theta_o = 0.006f;
const real tau_v1minus = 10.0f;
const real tau_v2minus = 1150.0f;
const real tau_vplus = 1.4506f;
const real tau_w1minus = 140.0f;
const real tau_w2minus = 6.25f;
const real k_wminus = 65.0f;
const real u_wminus = 0.015f;
const real tau_wplus = 326.0f;
const real tau_fi = 0.105f;
const real tau_o1 = 400.0f;
const real tau_o2 = 6.0f;
const real tau_so1 = 30.0181f;
const real tau_so2 = 0.9957f;
const real k_so = 2.0458f;
const real u_so = 0.65f;
const real tau_s1 = 2.7342f;
const real tau_s2 = 16.0f;
const real k_s = 2.0994f;
const real u_s = 0.9087f;
const real tau_si = 1.8875f;
const real tau_winf = 0.175f;
const real w_infstar = 0.9f;

#endif // PB

#ifdef TNNP

const real u_o = 0.0f;
const real u_u = 1.58f;
const real theta_v = 0.3f;
const real theta_w = 0.015f;
const real theta_vminus = 0.015f;
const real theta_o = 0.006f;
const real tau_v1minus = 60.0f;
const real tau_v2minus = 1150.0f;
const real tau_vplus = 1.4506f;
const real tau_w1minus = 70.0f;
const real tau_w2minus = 20.0f;
const real k_wminus = 65.0f;
const real u_wminus = 0.03f;
const real tau_wplus = 280.0f;
const real tau_fi = 0.11f;
const real tau_o1 = 6.0f;
const real tau_o2 = 6.0f;
const real tau_so1 = 43.0f;
const real tau_so2 = 0.2f;
const real k_so = 2.0f;
const real u_so = 0.65f;
const real tau_s1 = 2.7342f;
const real tau_s2 = 3.0f;
const real k_s = 2.0994f;
const real u_s = 0.9087f;
const real tau_si = 2.8723f;
const real tau_winf = 0.07f;
const real w_infstar = 0.94f;

#endif // TNNP

#endif // not USE_CUDA

void solveMonodomainMV(const SimulationConfig *config, Measurement *measurement, const real *time_array)
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