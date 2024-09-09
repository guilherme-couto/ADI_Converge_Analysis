#ifndef SIMULATION_CONFIG_H
#define SIMULATION_CONFIG_H

#include "include.h"

#ifdef SERIAL
const int L = 1;     // Length of each side (cm)
const real totalTime = 0.1f; // Total time (ms)

#ifdef MONODOMAIN
// Stimulation parameters
const real stimStrength = 100.0f;

// S1
const real stim1Begin = 0.0f;    // S1 start time -> ms
const real stim1Duration = 2.0f; // Stimulation duration -> ms
const real stim1xLimit = 0.2f;   // Stimulation x limit -> cm
const real stim1yLimit = 2.0f;   // Stimulation y limit -> cm ( = L)

// S2
const real stim2Begin = 120.0f;  // S2 start time -> ms
const real stim2Duration = 2.0f; // Stimulation duration -> ms
const real stim2xMax = 1.0f;     // Stimulation x max -> cm
const real stim2yMax = 1.0f;     // Stimulation y max -> cm
const real stim2xMin = 0.0f;     // Stimulation x min -> cm
const real stim2yMin = 0.0f;     // Stimulation y min -> cm

// Initial conditions
#ifdef AFHN
const real V_init = 0.0f;
const real W_init = 0.0f;
#endif // AFHN
#endif                           // MONODOMAIN
#endif                           // SERIAL

#ifdef GPU
const __constant__ int L = 5; // Length of each side (cm)
const real totalTime = 1000.0f;        // Total time (ms)

#ifdef SAVE_FRAMES
// Save frames parameters
const int frameSaveRate = 300;
#endif // SAVE_FRAMES

#ifdef MONODOMAIN
// Stimulation parameters
const __constant__ int numberOfStimuli = 2; // Number of stimuli
const real stimuliStrength = -38.0f;          // Stimulation strength -> (amplitude)
const real stimuliDuration = 2.0f;          // Stimulation duration -> ms
const real stimuliBegin[] = {0.0f, 250.0f}; // Stimuli begin time -> ms
const real stimulixMax[] = {0.2f, 2.5f};    // Stimuli x max -> cm
const real stimulixMin[] = {0.0f, 0.0f};    // Stimuli x min -> cm
const real stimuliyMax[] = {5.0f, 2.5f};    // Stimuli y max -> cm
const real stimuliyMin[] = {0.0f, 0.0f};    // Stimuli y min -> cm

// Initial conditions
#ifdef AFHN
const real V_init = 0.0f;
const real W_init = 0.0f;
#endif // AFHN

#ifdef TT2
#ifdef EPI
// Initial conditions for EPI cells
// from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3263775/
const real V_init = -85.23f;       // Initial membrane potential -> mV
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
const real V_init = -86.2f;       // Initial membrane potential -> mV
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
#endif                            // MONODOMAIN
#endif                            // GPU

#endif // SIMULATION_CONFIG_H