#ifndef INCLUDE_H
#define INCLUDE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <omp.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// ANSI color codes
#define RESET "\033[0m"
#define BLUE "\033[34m"
#define GREEN "\033[32m"
#define YELLOW "\033[33m"
#define RED "\033[31m"

// Macros for printing messages
#define INFOMSG(msg, ...) printf(BLUE "[i] " RESET msg, ##__VA_ARGS__)
#define SUCCESSMSG(msg, ...) printf(GREEN "[+] " RESET msg, ##__VA_ARGS__)
#define WARNINGMSG(msg, ...) printf(YELLOW "[!] " RESET msg, ##__VA_ARGS__)
#define ERRORMSG(msg, ...) printf(RED "[x] " RESET msg, ##__VA_ARGS__)
#define DEBUGMSG(msg) printf(YELLOW "[.] " RESET "Line number %d in file %s\n", __LINE__, __FILE__)

// Define maximum string size
#define MAX_STRING_SIZE 200

// Define block size for GPU
#define BLOCK_SIZE 256

// Convert CM to UM
#define CM_TO_UM(x) ((int)(x * 1.0e4))

// Define real type via compile command line (-D{OPTION}, USE_const real or USE_FLOAT)
#ifdef USE_DOUBLE
typedef double real;
#define REAL_TYPE "double"
#define STRTOREAL strtod
#define FSCANF_REAL "%le"
#else
typedef float real;
#define REAL_TYPE "float"
#define STRTOREAL strtof
#define FSCANF_REAL "%e"
#endif

#ifdef SERIAL
#define EXECUTION_TYPE "SERIAL"
#define RUNSIMULATION runSimulationSerial
#endif // SERIAL
#ifdef GPU
#define EXECUTION_TYPE "GPU"
#define RUNSIMULATION runSimulationGPU
#endif // GPU

// Define cell model via compile command line (-D{OPTION}, AFHN, TT2 or MV)
// if cell model has phenotype, define it as well
#ifdef AFHN
#define CELL_MODEL "AFHN"
#endif // AFHN

#ifdef TT2
#define CELL_MODEL "TT2"
// Options: ENDO, MCELL, EPI
#if !defined(MCELL) && !defined(EPI)
#define ENDO
#endif // not MCELL && not EPI
#endif // TT2

#ifdef MV
#define CELL_MODEL "MV"
// Options: ENDO, M, EPI, PB, TNNP
#if !defined(MCELL) && !defined(EPI) && !defined(PB) && !defined(TNNP)
#define ENDO
#endif // not MCELL && not EPI && not PB && not TNNP
#endif // MV

// Define CUDA error checking
#ifdef GPU
#define CUDA_CALL(call)                                                                                   \
    do                                                                                                    \
    {                                                                                                     \
        cudaError_t error = call;                                                                         \
        if (error != cudaSuccess)                                                                         \
        {                                                                                                 \
            fprintf(stderr, "CUDA error at %s:%d - %s\n", __FILE__, __LINE__, cudaGetErrorString(error)); \
            exit(EXIT_FAILURE);                                                                           \
        }                                                                                                 \
    } while (0)
#endif // GPU

// Define problem via compile command line (-D{OPTION}):
// forcing only for CONVERGENCE_ANALYSIS_FORCING_TERM
// LINMONO -> Adapted monodomain with linear reaction (2D)
//            chi*Cm*dv/dt = sigma*Lap(v) - chi*G*v + forcing
//            Boundaries: Neumann
//
// MONODOMAIN with AFHN -> Monodomain with adapted FitzHugh-Nagumo (2D)
//                         { dv/dt = (sigma/(chi*Cm))*Lap(v) - RHS_v/Cm [+ forcing/(chi*Cm)]
//                         { dw/dt = RHS_w
//                         RHS_v = (G*v*(1.0-(v/vth)) * (1.0-(v/vp))) + (eta1*v*w)
//                         RHS_w = eta2 * ((v/vp)-(eta3*w))
//                         Boundaries: Neumann
//
// CABLEEQ -> Cable equation with adapted FitzHugh-Nagumo (1D)
//            { dv/dt = (sigma/(chi*Cm))*d²v/dx² - RHS_v/Cm
//            { dw/dt = RHS_w
//            RHS_v = (G*v*(1.0-(v/vth)) * (1.0-(v/vp))) + (eta1*v*w)
//            RHS_w = eta2 * ((v/vp)-(eta3*w))
//            Boundaries: Neumann

#ifdef LINMONO
#define PROBLEM "LINMONO"
#endif // LINMONO
#ifdef MONODOMAIN
#define PROBLEM "MONODOMAIN"
#endif // MONODOMAIN
#ifdef CABLEEQ
#define PROBLEM "CABLEEQ"
#endif // CABLEEQ

// Define aux structures for MONODOMAIN
#if defined(MONODOMAIN) || defined(CABLEEQ)
typedef struct
{
    real strength;
    real begin;
    real duration;
    int xMaxDisc;
    int xMinDisc;
    int yMaxDisc;
    int yMinDisc;
} Stimulus;
#endif // MONODOMAIN || CABLEEQ

// Define method via compile command line (-D{OPTION}):
// ADI -> Alternating Direction Implicit
// OS-ADI -> Operator Splitting ADI
// SSI-ADI -> Second Order Semi Implicit ADI
// theta-ADI -> theta method with ADI
// theta-RK2 -> theta method with RK2 (only for CABLEEQ)
// FE -> Forward Euler

#ifdef ADI
#define METHOD "ADI"
#endif // ADI
#ifdef OSADI
#define METHOD "OS-ADI"
#endif // OSADI
#ifdef SSIADI
#define METHOD "SSI-ADI"
#endif // SSIADI
#ifdef THETASSIADI
#define METHOD "theta-SSI-ADI"
#endif // THETASSIADI
#ifdef THETASSIRK2
#define METHOD "theta-SSI-RK2"
#endif // THETARK2
#ifdef FE
#define METHOD "FE"
#endif // FE

// If defined SERIAL, constants are defined only as const for CPU
#ifdef SERIAL
const real _pi = 3.14159265358979323846f;

#ifdef LINMONO
const real G = 1.0f;     // omega^-1 * cm^-2
const real sigma = 1.0f; // omega^-1 * cm^-1
const real chi = 1.0f;   // cm^-1
const real Cm = 1.0f;    // mF * cm^-2
#endif                   // LINMONO
#ifdef DIFFREAC
const real sigma = 1.0f;
#endif // DIFFREAC
#ifdef DIFF
const real sigma = 1.0f;
#endif // DIFF
#if defined(MONODOMAIN) || defined(CABLEEQ)
#if defined(CONVERGENCE_ANALYSIS_FORCING_TERM) && defined(AFHN)
const real sigma = 1.0f; // omega^-1 * cm^-1
const real chi = 1.0f;   // cm^-1
const real Cm = 1.0f;    // mF * cm^-2

const real G = 1.0f;    // omega^-1 * cm^-2
const real eta1 = 1.0f; // omega^-1 * cm^-1
const real eta2 = 1.0f; // dimensionless
const real eta3 = 1.0f; // dimensionless
const real vth = 1.0f;  // mV
const real vp = 1.0f;   // mV
#elif defined(AFHN)
// Model parameters - Based on Gerardo_Giorda 2007
const real sigma = 1.2e-3f; // omega^-1 * cm^-1
const real chi = 1.0e3f;    // cm^-1
const real Cm = 1.0e-3f;    // mF * cm^-2

const real G = 1.5f;      // omega^-1 * cm^-2
const real eta1 = 4.4f;   // omega^-1 * cm^-1
const real eta2 = 0.012f; // dimensionless
const real eta3 = 1.0f;   // dimensionless
const real vth = 13.0f;   // mV
const real vp = 100.0f;   // mV
#endif // AFHN

#ifdef TT2
// Model parameters - Based on Ten Tusscher 2006 (https://journals.physiology.org/doi/full/10.1152/ajpheart.00109.2006)
// from https://tbb.bio.uu.nl/khwjtuss/SourceCodes/HVM2/Source/Main.cc - ten Tusscher code
// https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3263775/ - Benchmark
// and https://github.com/rsachetto/MonoAlg3D_C/blob/master/src/models_library/ten_tusscher/ten_tusscher_2006_RS_CPU.c - Sachetto MonoAlg3D
// --------------------------------------------------------------------------------------------------------------------------------------------------------*/
const real chi = 1400.0f;  // Surface-to-volume ratio (cm^-1)
const real Cm = 0.185f;    // Membrane capacitance (uF/cm^2)
const real sigma = 1.171f; // Diffusion coefficient (cm^2/s)

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
#endif                    // EPI || MCELL
#ifdef ENDO
const real G_to = 0.073f; // Maximal I_to (transient outward potassium current) conductance -> nS/pF (endo cells)
#endif                    // ENDO
const real G_Kr = 0.153f; // Maximal I_Kr (rapidly activating delayed rectifier potassium current) conductance -> nS/pF
#if defined(EPI) || defined(ENDO)
const real G_Ks = 0.392f; // Maximal I_Ks (slowly activating delayed rectifier potassium current) conductance -> nS/pF (epi and endo cells)
#endif                    // EPI || ENDO
#ifdef MCELL
const real G_Ks = 0.098f;        // Maximal I_Ks (slowly activating delayed rectifier potassium current) conductance -> nS/pF (M cells)
#endif                           // MCELL
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
#endif                         // TT2

#ifdef MV
// Model definition https://www.sciencedirect.com/science/article/pii/S0022519308001690?via%3Dihub
const real Dtilde = 1.171f;
const real chi = 1400.0f;

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
#endif // MV
#endif // MONODOMAIN || CABLEEQ
#endif // SERIAL

// If defined GPU, constants are defined as const for CPU and __constant__ for GPU
#ifdef GPU
const __constant__ real _pi = 3.14159265358979323846f;

#ifdef LINMONO
const __constant__ real G = 1.0f;     // omega^-1 * cm^-2
const __constant__ real sigma = 1.0f; // omega^-1 * cm^-1
const __constant__ real chi = 1.0f;   // cm^-1
const __constant__ real Cm = 1.0f;    // mF * cm^-2
#endif                                // LINMONO
#ifdef DIFFREAC
const __constant__ real sigma = 1.0f;
#endif // DIFFREAC
#ifdef DIFF
const __constant__ real sigma = 1.0f;
#endif // DIFF
#ifdef MONODOMAIN
#if defined(CONVERGENCE_ANALYSIS_FORCING_TERM) && defined(AFHN)
const __constant__ real sigma = 1.0f; // omega^-1 * cm^-1
const __constant__ real chi = 1.0f;   // cm^-1
const __constant__ real Cm = 1.0f;    // mF * cm^-2

const __constant__ real G = 1.0f;    // omega^-1 * cm^-2
const __constant__ real eta1 = 1.0f; // omega^-1 * cm^-1
const __constant__ real eta2 = 1.0f; // dimensionless
const __constant__ real eta3 = 1.0f; // dimensionless
const __constant__ real vth = 1.0f;  // mV
const __constant__ real vp = 1.0f;   // mV
#elif defined(AFHN)
// Model parameters - Based on Gerardo_Giorda 2007
const __constant__ real sigma = 1.2e-3f; // omega^-1 * cm^-1
const __constant__ real chi = 1.0e3f;    // cm^-1
const __constant__ real Cm = 1.0e-3f;    // mF * cm^-2

const __constant__ real G = 1.5f;      // omega^-1 * cm^-2
const __constant__ real eta1 = 4.4f;   // omega^-1 * cm^-1
const __constant__ real eta2 = 0.012f; // dimensionless
const __constant__ real eta3 = 1.0f;   // dimensionless
const __constant__ real vth = 13.0f;   // mV
const __constant__ real vp = 100.0f;   // mV
#endif // AFHN

#ifdef TT2
// Model parameters - Based on Ten Tusscher 2006 (https://journals.physiology.org/doi/full/10.1152/ajpheart.00109.2006)
// from https://tbb.bio.uu.nl/khwjtuss/SourceCodes/HVM2/Source/Main.cc - ten Tusscher code
// https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3263775/ - Benchmark
// and https://github.com/rsachetto/MonoAlg3D_C/blob/master/src/models_library/ten_tusscher/ten_tusscher_2006_RS_CPU.c - Sachetto MonoAlg3D
// --------------------------------------------------------------------------------------------------------------------------------------------------------*/
const __constant__ real chi = 1400.0f;  // Surface-to-volume ratio (cm^-1)
const __constant__ real Cm = 0.185f;    // Membrane capacitance (uF/cm^2)
const __constant__ real sigma = 1.171f; // Diffusion coefficient (cm^2/s)

const __constant__ real R = 8314.472f;      // Universal gas constant (J/(kmol*K))
const __constant__ real T = 310.0f;         // Temperature (K)
const __constant__ real F = 96485.3415f;    // Faraday's constant (C/mol)
const __constant__ real RTONF = 26.713761f; // R*T/F
const __constant__ real FONRT = 0.037434f;  // F/(R*T)

// Intracellular volumes
const __constant__ real V_C = 0.016404f;    // Cellular volume -> (???) [16404 um^3]
const __constant__ real V_SR = 0.001094f;   // Sarcoplasmic reticulum volume -> (???) [1094 um^3]
const __constant__ real V_SS = 0.00005468f; // Subsarcolemmal space volume -> (???) [54.68 um^3]

// External concentrations
const __constant__ real K_o = 5.4f;    // Extracellular potassium (K+) concentration -> mM
const __constant__ real Na_o = 140.0f; // Extracellular sodium (Na+) concentration -> mM
const __constant__ real Ca_o = 2.0f;   // Extracellular calcium (Ca++) concentration -> mM

// Parameters for currents
const __constant__ real G_Na = 14.838f; // Maximal I_Na (sodium current) conductance -> nS/pF
const __constant__ real G_K1 = 5.405f;  // Maximal I_K1 (late rectifier potassium current) conductance -> nS/pF
#if defined(EPI) || defined(MCELL)
const __constant__ real G_to = 0.294f; // Maximal I_to (transient outward potassium current) conductance -> nS/pF (epi and M cells)
#endif                                 // EPI || MCELL
#ifdef ENDO
const __constant__ real G_to = 0.073f; // Maximal I_to (transient outward potassium current) conductance -> nS/pF (endo cells)
#endif                                 // ENDO
const __constant__ real G_Kr = 0.153f; // Maximal I_Kr (rapidly activating delayed rectifier potassium current) conductance -> nS/pF
#if defined(EPI) || defined(ENDO)
const __constant__ real G_Ks = 0.392f; // Maximal I_Ks (slowly activating delayed rectifier potassium current) conductance -> nS/pF (epi and endo cells)
#endif                                 // EPI || ENDO
#ifdef MCELL
const __constant__ real G_Ks = 0.098f;        // Maximal I_Ks (slowly activating delayed rectifier potassium current) conductance -> nS/pF (M cells)
#endif                                        // MCELL
const __constant__ real p_KNa = 0.03f;        // Relative I_Ks permeability to Na+ over K+ -> dimensionless
const __constant__ real G_CaL = 3.98e-5f;     // Maximal I_CaL (L-type calcium current) conductance -> cm/ms/uF
const __constant__ real k_NaCa = 1000.0f;     // Maximal I_NaCa (Na+/Ca++ exchanger current) -> pA/pF
const __constant__ real gamma_I_NaCa = 0.35f; // Voltage dependence parameter of I_NaCa -> dimensionless
const __constant__ real K_mCa = 1.38f;        // Half-saturation constant of I_NaCa for intracellular Ca++ -> mM
const __constant__ real K_mNa_i = 87.5f;      // Half-saturation constant of I_NaCa for intracellular Na+ -> mM
const __constant__ real k_sat = 0.1f;         // Saturation factor for I_NaCa -> dimensionless
const __constant__ real alpha = 2.5f;         // Factor enhancing outward nature of I_NaCa -> dimensionless
const __constant__ real P_NaK = 2.724f;       // Maximal I_NaK (Na+/K+ pump current) -> pA/pF
const __constant__ real K_mK = 1.0f;          // Half-saturation constant of I_NaK for Ko -> mM
const __constant__ real K_mNa = 40.0f;        // Half-saturation constant of I_NaK for intracellular Na+ -> mM
const __constant__ real G_pK = 0.0146f;       // Maximal I_pK (plateau potassium current) conductance -> nS/pF
const __constant__ real G_pCa = 0.1238f;      // Maximal I_pCa (plateau calcium current) conductance -> nS/pF
const __constant__ real K_pCa = 0.0005f;      // Half-saturation constant of I_pCa for intracellular Ca++ -> mM
const __constant__ real G_bNa = 0.00029f;     // Maximal I_bNa (sodium background current) conductance -> nS/pF
const __constant__ real G_bCa = 0.000592f;    // Maximal I_bCa (calcium background current) conductance -> nS/pF

// Intracellular calcium flux dynamics
const __constant__ real V_maxup = 0.006375f; // Maximal I_up -> mM/ms
const __constant__ real K_up = 0.00025f;     // Half-saturation constant of I_up -> mM
const __constant__ real V_rel = 0.102f;      // Maximal I_rel conductance -> mM/ms
const __constant__ real k1_prime = 0.15f;    // R to O and RI to I I_rel transition rate -> mM^-2*ms^-1
const __constant__ real k2_prime = 0.045f;   // O to I  and R to RI I_rel transition rate -> mM^-1*ms^-1
const __constant__ real k3 = 0.06f;          // O to R and I to RI I_rel transition rate -> ms^-1
const __constant__ real k4 = 0.005f;         // I to O and RI to I I_rel transition rate -> ms^-1
const __constant__ real EC = 1.5f;           // Half-saturation constant of k_Ca_SR -> mM
const __constant__ real max_SR = 2.5f;       // Maximum value of k_Ca_SR -> dimensionless
const __constant__ real min_SR = 1.0f;       // Minimum value of k_Ca_SR -> dimensionless
const __constant__ real V_leak = 0.00036f;   // Maximal I_leak conductance -> mM/ms
const __constant__ real V_xfer = 0.0038f;    // Maximal I_xfer conductance -> mM/ms

// Calcium buffering dynamics
const __constant__ real bufC = 0.2f;        // Total cytoplasmic buffer concentration -> mM
const __constant__ real K_bufC = 0.001f;    // Half-saturation constant of cytoplasmic buffers -> mM
const __constant__ real bufSR = 10.0f;      // Total sarcoplasmic reticulum buffer concentration -> mM
const __constant__ real K_bufSR = 0.3f;     // Half-saturation constant of sarcoplasmic reticulum buffers -> mM
const __constant__ real bufSS = 0.4f;       // Total subspace buffer concentration -> mM
const __constant__ real K_bufSS = 0.00025f; // Half-saturation constant of subspace buffer -> mM
#endif                                      // TT2

#ifdef MV
// Model definition https://www.sciencedirect.com/science/article/pii/S0022519308001690?via%3Dihub
const __constant__ real Dtilde = 1.171f;
const __constant__ real chi = 1400.0f;

#ifdef EPI
const __constant__ real u_o = 0.0f;
const __constant__ real u_u = 1.55f;
const __constant__ real theta_v = 0.3f;
const __constant__ real theta_w = 0.13f;
const __constant__ real theta_vminus = 0.006f;
const __constant__ real theta_o = 0.006f;
const __constant__ real tau_v1minus = 60.0f;
const __constant__ real tau_v2minus = 1150.0f;
const __constant__ real tau_vplus = 1.4506f;
const __constant__ real tau_w1minus = 60.0f;
const __constant__ real tau_w2minus = 15.0f;
const __constant__ real k_wminus = 65.0f;
const __constant__ real u_wminus = 0.03f;
const __constant__ real tau_wplus = 200.0f;
const __constant__ real tau_fi = 0.11f;
const __constant__ real tau_o1 = 400.0f;
const __constant__ real tau_o2 = 6.0f;
const __constant__ real tau_so1 = 30.0181f;
const __constant__ real tau_so2 = 0.9957f;
const __constant__ real k_so = 2.0458f;
const __constant__ real u_so = 0.65f;
const __constant__ real tau_s1 = 2.7342f;
const __constant__ real tau_s2 = 16.0f;
const __constant__ real k_s = 2.0994f;
const __constant__ real u_s = 0.9087f;
const __constant__ real tau_si = 1.8875f;
const __constant__ real tau_winf = 0.07f;
const __constant__ real w_infstar = 0.94f;
#endif // EPI
#ifdef ENDO
const __constant__ real u_o = 0.0f;
const __constant__ real u_u = 1.56f;
const __constant__ real theta_v = 0.3f;
const __constant__ real theta_w = 0.13f;
const __constant__ real theta_vminus = 0.2f;
const __constant__ real theta_o = 0.006f;
const __constant__ real tau_v1minus = 75.0f;
const __constant__ real tau_v2minus = 10.0f;
const __constant__ real tau_vplus = 1.4506f;
const __constant__ real tau_w1minus = 6.0f;
const __constant__ real tau_w2minus = 140.0f;
const __constant__ real k_wminus = 200.0f;
const __constant__ real u_wminus = 0.016f;
const __constant__ real tau_wplus = 280.0f;
const __constant__ real tau_fi = 0.1f;
const __constant__ real tau_o1 = 470.0f;
const __constant__ real tau_o2 = 6.0f;
const __constant__ real tau_so1 = 40.0f;
const __constant__ real tau_so2 = 1.2f;
const __constant__ real k_so = 2.0f;
const __constant__ real u_so = 0.65f;
const __constant__ real tau_s1 = 2.7342f;
const __constant__ real tau_s2 = 2.0f;
const __constant__ real k_s = 2.0994f;
const __constant__ real u_s = 0.9087f;
const __constant__ real tau_si = 2.9013f;
const __constant__ real tau_winf = 0.0273f;
const __constant__ real w_infstar = 0.78f;
#endif // ENDO
#ifdef MCELL
const __constant__ real u_o = 0.0f;
const __constant__ real u_u = 1.61f;
const __constant__ real theta_v = 0.3f;
const __constant__ real theta_w = 0.13f;
const __constant__ real theta_vminus = 0.1f;
const __constant__ real theta_o = 0.005f;
const __constant__ real tau_v1minus = 80.0f;
const __constant__ real tau_v2minus = 1.4506f;
const __constant__ real tau_vplus = 1.4506f;
const __constant__ real tau_w1minus = 70.0f;
const __constant__ real tau_w2minus = 8.0f;
const __constant__ real k_wminus = 200.0f;
const __constant__ real u_wminus = 0.016f;
const __constant__ real tau_wplus = 280.0f;
const __constant__ real tau_fi = 0.078f;
const __constant__ real tau_o1 = 410.0f;
const __constant__ real tau_o2 = 7.0f;
const __constant__ real tau_so1 = 91.0f;
const __constant__ real tau_so2 = 0.8f;
const __constant__ real k_so = 2.1f;
const __constant__ real u_so = 0.6f;
const __constant__ real tau_s1 = 2.7342f;
const __constant__ real tau_s2 = 4.0f;
const __constant__ real k_s = 2.0994f;
const __constant__ real u_s = 0.9087f;
const __constant__ real tau_si = 3.3849f;
const __constant__ real tau_winf = 0.01f;
const __constant__ real w_infstar = 0.5f;
#endif // MCELL
#ifdef PB
const __constant__ real u_o = 0.0f;
const __constant__ real u_u = 1.45f;
const __constant__ real theta_v = 0.35f;
const __constant__ real theta_w = 0.13f;
const __constant__ real theta_vminus = 0.175f;
const __constant__ real theta_o = 0.006f;
const __constant__ real tau_v1minus = 10.0f;
const __constant__ real tau_v2minus = 1150.0f;
const __constant__ real tau_vplus = 1.4506f;
const __constant__ real tau_w1minus = 140.0f;
const __constant__ real tau_w2minus = 6.25f;
const __constant__ real k_wminus = 65.0f;
const __constant__ real u_wminus = 0.015f;
const __constant__ real tau_wplus = 326.0f;
const __constant__ real tau_fi = 0.105f;
const __constant__ real tau_o1 = 400.0f;
const __constant__ real tau_o2 = 6.0f;
const __constant__ real tau_so1 = 30.0181f;
const __constant__ real tau_so2 = 0.9957f;
const __constant__ real k_so = 2.0458f;
const __constant__ real u_so = 0.65f;
const __constant__ real tau_s1 = 2.7342f;
const __constant__ real tau_s2 = 16.0f;
const __constant__ real k_s = 2.0994f;
const __constant__ real u_s = 0.9087f;
const __constant__ real tau_si = 1.8875f;
const __constant__ real tau_winf = 0.175f;
const __constant__ real w_infstar = 0.9f;
#endif // PB
#ifdef TNNP
const __constant__ real u_o = 0.0f;
const __constant__ real u_u = 1.58f;
const __constant__ real theta_v = 0.3f;
const __constant__ real theta_w = 0.015f;
const __constant__ real theta_vminus = 0.015f;
const __constant__ real theta_o = 0.006f;
const __constant__ real tau_v1minus = 60.0f;
const __constant__ real tau_v2minus = 1150.0f;
const __constant__ real tau_vplus = 1.4506f;
const __constant__ real tau_w1minus = 70.0f;
const __constant__ real tau_w2minus = 20.0f;
const __constant__ real k_wminus = 65.0f;
const __constant__ real u_wminus = 0.03f;
const __constant__ real tau_wplus = 280.0f;
const __constant__ real tau_fi = 0.11f;
const __constant__ real tau_o1 = 6.0f;
const __constant__ real tau_o2 = 6.0f;
const __constant__ real tau_so1 = 43.0f;
const __constant__ real tau_so2 = 0.2f;
const __constant__ real k_so = 2.0f;
const __constant__ real u_so = 0.65f;
const __constant__ real tau_s1 = 2.7342f;
const __constant__ real tau_s2 = 3.0f;
const __constant__ real k_s = 2.0994f;
const __constant__ real u_s = 0.9087f;
const __constant__ real tau_si = 2.8723f;
const __constant__ real tau_winf = 0.07f;
const __constant__ real w_infstar = 0.94f;
#endif // TNNP
#endif // MV
#endif // MONODOMAIN
#endif // GPU

#endif // INCLUDE_H