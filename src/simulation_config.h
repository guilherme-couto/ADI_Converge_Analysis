#ifndef SIMULATION_CONFIG_H
#define SIMULATION_CONFIG_H

#include "include.h"

#ifdef SAVE_FRAMES
// const int frameSaveRate = 20000; // dt = 0.0001 ms saving each 2 ms
// const int frameSaveRate = 500; // dt = 0.001 ms saving each 2 ms
const int frameSaveRate = 200; // dt = 0.01 ms saving each 2 ms
#endif // SAVE_FRAMES

#ifndef CONVERGENCE_ANALYSIS_FORCING_TERM

#ifdef CABLEEQ

const real totalTime = 50.0f;      // Total time (ms)
const __constant__ real Lx = 5.0f; // Length in x (cm)

#else // if not CABLEEQ

const real totalTime = 100.0f;       // Total time (ms)
const __constant__ real Lx = 6.0f;  // Length in x (cm)
const __constant__ real Ly = 6.0f; // Length in y (cm)

#endif // CABLEEQ

#else // if CONVERGENCE_ANALYSIS_FORCING_TERM

const real totalTime = 0.1f;       // Total time (ms)
const __constant__ real Lx = 1.0f; // Length in x (cm)
const __constant__ real Ly = 1.0f; // Length in y (cm)

#endif // not CONVERGENCE_ANALYSIS_FORCING_TERM

#if defined(MONODOMAIN) || defined(CABLEEQ)

// Stimulation parameters
const __constant__ int numberOfStimuli = 2; // Number of stimuli

#ifdef AFHN
const real stimuliAmplitude = 100.0f; // Stimulation amplitude -> (amplitude)
#endif                               // AFHN

#ifdef TT2
const real stimuliAmplitude = 38.0f; // Stimulation amplitude -> (amplitude)
#endif                              // TT2

#ifdef MV
const real stimuliAmplitude = 1.0f; // Stimulation amplitude -> (amplitude)
#endif                             // MV

const real stimuliDuration = 2.0f;          // Stimulation duration -> ms
const real stimuliBegin[] = {0.0f, 340.0f}; // Stimuli begin time -> ms
const real stimulixMax[] = {0.2f, 3.0f};    // Stimuli x max -> cm
const real stimulixMin[] = {0.0f, 0.0f};    // Stimuli x min -> cm
const real stimuliyMax[] = {6.0f, 3.0f};    // Stimuli y max -> cm
const real stimuliyMin[] = {0.0f, 0.0f};    // Stimuli y min -> cm

// Initial conditions
#ifdef AFHN
const real Vm_init = 0.0f;
const real W_init = 0.0f;
#endif // AFHN

#ifdef TT2
#ifdef EPI
// Initial conditions for EPI cells
// from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3263775/
const real Vm_init = -85.23f;      // Initial membrane potential -> mV
const real X_r1_init = 0.00621f;   // Initial rapid time-dependent potassium current Xr1 gate -> dimensionless
const real X_r2_init = 0.4712f;    // Initial rapid time-dependent potassium current Xr2 gate -> dimensionless
const real X_s_init = 0.0095f;     // Initial slow time-dependent potassium current Xs gate -> dimensionless
const real m_init = 0.00172f;      // Initial fast sodium current m gate -> dimensionless
const real h_init = 0.7444f;       // Initial fast sodium current h gate -> dimensionless
const real j_init = 0.7045f;       // Initial fast sodium current j gate -> dimensionless
const real d_init = 3.373e-5f;     // Initial L-type calcium current d gate -> dimensionless
const real f_init = 0.7888f;       // Initial L-type calcium current f gate -> dimensionless
const real f2_init = 0.9755f;      // Initial L-type calcium current f2 gate -> dimensionless
const real fCaSS_init = 0.9953f;   // Initial L-type calcium current fCaSS gate -> dimensionless
const real s_init = 0.999998f;     // Initial transient outward current s gate -> dimensionless
const real r_init = 2.42e-8f;      // Initial transient outward current r gate -> dimensionless
const real Ca_i_init = 0.000126f;  // Initial intracellular Ca++ concentration -> mM
const real Ca_SR_init = 3.64f;     // Initial sarcoplasmic reticulum Ca++ concentration -> mM
const real Ca_SS_init = 0.00036f;  // Initial subspace Ca++ concentration -> mM
const real R_prime_init = 0.9073f; // Initial ryanodine receptor -> dimensionless
const real Na_i_init = 8.604f;     // Initial intracellular Na+ concentration -> mM
const real K_i_init = 136.89f;     // Initial intracellular K+ concentration -> mM
#endif                             // EPI
#if defined(ENDO) || defined(MCELL)
// Initial conditions for ENDO and MCELL cells
// from https://tbb.bio.uu.nl/khwjtuss/SourceCodes/HVM2/
const real Vm_init = -86.2f;      // Initial membrane potential -> mV
const real X_r1_init = 0.0f;      // Initial rapid time-dependent potassium current Xr1 gate -> dimensionless
const real X_r2_init = 1.0f;      // Initial rapid time-dependent potassium current Xr2 gate -> dimensionless
const real X_s_init = 0.0f;       // Initial slow time-dependent potassium current Xs gate -> dimensionless
const real m_init = 0.0f;         // Initial fast sodium current m gate -> dimensionless
const real h_init = 0.75f;        // Initial fast sodium current h gate -> dimensionless
const real j_init = 0.75f;        // Initial fast sodium current j gate -> dimensionless
const real d_init = 0.0f;         // Initial L-type calcium current d gate -> dimensionless
const real f_init = 1.0f;         // Initial L-type calcium current f gate -> dimensionless
const real f2_init = 1.0f;        // Initial L-type calcium current f2 gate -> dimensionless
const real fCaSS_init = 1.0f;     // Initial L-type calcium current fCaSS gate -> dimensionless
const real s_init = 1.0f;         // Initial transient outward current s gate -> dimensionless
const real r_init = 0.0f;         // Initial transient outward current r gate -> dimensionless
const real Ca_i_init = 0.00007f;  // Initial intracellular Ca++ concentration -> mM
const real Ca_SR_init = 1.3f;     // Initial sarcoplasmic reticulum Ca++ concentration -> mM
const real Ca_SS_init = 0.00007f; // Initial subspace Ca++ concentration -> mM
const real R_prime_init = 1.0f;   // Initial ryanodine receptor -> dimensionless
const real Na_i_init = 7.67f;     // Initial intracellular Na+ concentration -> mM
const real K_i_init = 138.3f;     // Initial intracellular K+ concentration -> mM
#endif                            // ENDO || MCELL
#endif                            // TT2

#ifdef MV
const real u_init = 0.0f;
const real v_init = 1.0f;
const real w_init = 1.0f;
const real s_init = 0.0f;
#endif // MV
#endif // MONODOMAIN || CABLEEQ

#endif // SIMULATION_CONFIG_H
