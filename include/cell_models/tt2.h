#ifndef TT2_MODEL_H
#define TT2_MODEL_H

#include "../core_definitions.h"

#ifdef TT2
#define CELL_MODEL "TT2"

// Options: ENDO, MCELL, EPI
#if !defined(MCELL) && !defined(EPI)
#define ENDO
#endif

// Model parameters - Based on Ten Tusscher 2006 (https://journals.physiology.org/doi/full/10.1152/ajpheart.00109.2006)
// from https://tbb.bio.uu.nl/khwjtuss/SourceCodes/HVM2/Source/Main.cc - ten Tusscher code
// https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3263775/ - Benchmark
// and https://github.com/rsachetto/MonoAlg3D_C/blob/master/src/models_library/ten_tusscher/ten_tusscher_2006_RS_CPU.c - Sachetto MonoAlg3D

#if defined(SERIAL) || defined(OPENMP)

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

#endif // SERIAL || OPENMP

#ifdef GPU

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

#endif // EPI || MCELL

#ifdef ENDO

const __constant__ real G_to = 0.073f; // Maximal I_to (transient outward potassium current) conductance -> nS/pF (endo cells)

#endif // ENDO

const __constant__ real G_Kr = 0.153f; // Maximal I_Kr (rapidly activating delayed rectifier potassium current) conductance -> nS/pF

#if defined(EPI) || defined(ENDO)

const __constant__ real G_Ks = 0.392f; // Maximal I_Ks (slowly activating delayed rectifier potassium current) conductance -> nS/pF (epi and endo cells)

#endif // EPI || ENDO

#ifdef MCELL

const __constant__ real G_Ks = 0.098f; // Maximal I_Ks (slowly activating delayed rectifier potassium current) conductance -> nS/pF (M cells)

#endif // MCELL

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

#endif // GPU

#endif // TT2

#endif // TT2_MODEL_H