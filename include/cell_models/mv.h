#ifndef MV_MODEL_H
#define MV_MODEL_H

#include "../core_definitions.h"

#ifdef MV
#define CELL_MODEL "MV"

// Options: ENDO, M, EPI, PB, TNNP
#if !defined(MCELL) && !defined(EPI) && !defined(PB) && !defined(TNNP)
#define ENDO
#endif

// Parameters - Based on MV model
// Model definition https://www.sciencedirect.com/science/article/pii/S0022519308001690?via%3Dihub

#if defined(SERIAL) || defined(OPENMP)

const real Dtilde = 0.65f * 1e-3; // cm^2/s

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

#endif // SERIAL || OPENMP

#ifdef GPU

const __constant__ real Dtilde = 0.65f * 1e-3; // cm^2/s

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

#endif // GPU

#endif // MV

#endif // MV_MODEL_H