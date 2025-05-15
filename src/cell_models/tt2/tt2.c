#include "tt2.h"

#ifndef USE_CUDA

// Model parameters - Based on Ten Tusscher 2006 (https://journals.physiology.org/doi/full/10.1152/ajpheart.00109.2006)
// from https://tbb.bio.uu.nl/khwjtuss/SourceCodes/HVM2/Source/Main.cc - ten Tusscher code
// https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3263775/ - Benchmark
// and https://github.com/rsachetto/MonoAlg3D_C/blob/master/src/models_library/ten_tusscher/ten_tusscher_2006_RS_CPU.c - Sachetto MonoAlg3D


// const real chi = 1400.0f;  // Surface-to-volume ratio (cm^-1)
// const real Cm = 0.185f;    // Membrane capacitance (uF/cm^2)
// const real sigma = 1.171f; // Diffusion coefficient (cm^2/s)

const real R = 8314.472f;      // Universal gas constant (J/(kmol*K))
const real T = 310.0f;         // Temperature (K)
const real F = 96485.3415f;    // Faraday's constant (C/mol)
const real RTONF = 26.713761f; // R*T/F
const real FONRT = 0.037434f;  // F/(R*T)

// Intracellular volumes
const real V_C = 0.016404f;    // Cellular volume -> (???) [16404 um^3]
const real V_SR = 0.001094f;   // Sarcoplasmic reticulum volume -> (???) [1094 um^3]
const real V_SS = 0.00005468f; // Subsarcolemmal space volume -> (???) [54.68 um^3]

// External concentrations
const real K_o = 5.4f;    // Extracellular potassium (K+) concentration -> mM
const real Na_o = 140.0f; // Extracellular sodium (Na+) concentration -> mM
const real Ca_o = 2.0f;   // Extracellular calcium (Ca++) concentration -> mM

// Parameters for currents
const real G_Na = 14.838f; // Maximal I_Na (sodium current) conductance -> nS/pF
const real G_K1 = 5.405f;  // Maximal I_K1 (late rectifier potassium current) conductance -> nS/pF

#if defined(EPI) || defined(MCELL)

const real G_to = 0.294f; // Maximal I_to (transient outward potassium current) conductance -> nS/pF (epi and M cells)

#endif // EPI || MCELL

#ifdef ENDO

const real G_to = 0.073f; // Maximal I_to (transient outward potassium current) conductance -> nS/pF (endo cells)

#endif // ENDO

const real G_Kr = 0.153f; // Maximal I_Kr (rapidly activating delayed rectifier potassium current) conductance -> nS/pF

#if defined(EPI) || defined(ENDO)

const real G_Ks = 0.392f; // Maximal I_Ks (slowly activating delayed rectifier potassium current) conductance -> nS/pF (epi and endo cells)

#endif // EPI || ENDO

#ifdef MCELL

const real G_Ks = 0.098f; // Maximal I_Ks (slowly activating delayed rectifier potassium current) conductance -> nS/pF (M cells)

#endif // MCELL

const real p_KNa = 0.03f;        // Relative I_Ks permeability to Na+ over K+ -> dimensionless
const real G_CaL = 3.98e-5f;     // Maximal I_CaL (L-type calcium current) conductance -> cm/ms/uF
const real k_NaCa = 1000.0f;     // Maximal I_NaCa (Na+/Ca++ exchanger current) -> pA/pF
const real gamma_I_NaCa = 0.35f; // Voltage dependence parameter of I_NaCa -> dimensionless
const real K_mCa = 1.38f;        // Half-saturation constant of I_NaCa for intracellular Ca++ -> mM
const real K_mNa_i = 87.5f;      // Half-saturation constant of I_NaCa for intracellular Na+ -> mM
const real k_sat = 0.1f;         // Saturation factor for I_NaCa -> dimensionless
const real alpha = 2.5f;         // Factor enhancing outward nature of I_NaCa -> dimensionless
const real P_NaK = 2.724f;       // Maximal I_NaK (Na+/K+ pump current) -> pA/pF
const real K_mK = 1.0f;          // Half-saturation constant of I_NaK for Ko -> mM
const real K_mNa = 40.0f;        // Half-saturation constant of I_NaK for intracellular Na+ -> mM
const real G_pK = 0.0146f;       // Maximal I_pK (plateau potassium current) conductance -> nS/pF
const real G_pCa = 0.1238f;      // Maximal I_pCa (plateau calcium current) conductance -> nS/pF
const real K_pCa = 0.0005f;      // Half-saturation constant of I_pCa for intracellular Ca++ -> mM
const real G_bNa = 0.00029f;     // Maximal I_bNa (sodium background current) conductance -> nS/pF
const real G_bCa = 0.000592f;    // Maximal I_bCa (calcium background current) conductance -> nS/pF

// Intracellular calcium flux dynamics
const real V_maxup = 0.006375f; // Maximal I_up -> mM/ms
const real K_up = 0.00025f;     // Half-saturation constant of I_up -> mM
const real V_rel = 0.102f;      // Maximal I_rel conductance -> mM/ms
const real k1_prime = 0.15f;    // R to O and RI to I I_rel transition rate -> mM^-2*ms^-1
const real k2_prime = 0.045f;   // O to I  and R to RI I_rel transition rate -> mM^-1*ms^-1
const real k3 = 0.06f;          // O to R and I to RI I_rel transition rate -> ms^-1
const real k4 = 0.005f;         // I to O and RI to I I_rel transition rate -> ms^-1
const real EC = 1.5f;           // Half-saturation constant of k_Ca_SR -> mM
const real max_SR = 2.5f;       // Maximum value of k_Ca_SR -> dimensionless
const real min_SR = 1.0f;       // Minimum value of k_Ca_SR -> dimensionless
const real V_leak = 0.00036f;   // Maximal I_leak conductance -> mM/ms
const real V_xfer = 0.0038f;    // Maximal I_xfer conductance -> mM/ms

// Calcium buffering dynamics
const real bufC = 0.2f;        // Total cytoplasmic buffer concentration -> mM
const real K_bufC = 0.001f;    // Half-saturation constant of cytoplasmic buffers -> mM
const real bufSR = 10.0f;      // Total sarcoplasmic reticulum buffer concentration -> mM
const real K_bufSR = 0.3f;     // Half-saturation constant of sarcoplasmic reticulum buffers -> mM
const real bufSS = 0.4f;       // Total subspace buffer concentration -> mM
const real K_bufSS = 0.00025f; // Half-saturation constant of subspace buffer -> mM

#endif // not USE_CUDA

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