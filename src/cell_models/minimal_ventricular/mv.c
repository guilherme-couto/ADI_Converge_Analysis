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

// MV solver functions
static inline void initialize_MV(real *restrict Vm, real *restrict sV, const int Nx, const int Ny)
{
    // Initial conditions
    const real u_init = 0.0f;
    const real v_init = 1.0f;
    const real w_init = 1.0f;
    const real s_init = 0.0f;

    // Initialize Vm and W with the initial condition
    for (int idx = 0, idx_sv = 0; idx < Nx * Ny; idx++, idx_sv = idx * MV_NSV)
    {
        Vm[idx] = u_init;
        sV[idx_sv] = v_init;
        sV[idx_sv + 1] = w_init;
        sV[idx_sv + 2] = s_init;
    }
}

static inline void get_actual_sV_MV(real *restrict actualsV, const real *restrict sV, const int idx)
{
    int idx_sv = idx * MV_NSV;
    actualsV[0] = sV[idx_sv];       // v
    actualsV[1] = sV[idx_sv + 1];   // w
    actualsV[2] = sV[idx_sv + 2];   // s
}

static inline const real compute_dVmdt_MV(const real Vm, const real *restrict sV)
{
    // Extract state variables
    const real v = sV[0];
    const real w = sV[1];
    const real s = sV[2];

    // Auxiliary variables
    const real Htheta_w = (Vm - theta_w > 0.0f) ? 1.0f : 0.0f;
    const real Htheta_o = (Vm - theta_o > 0.0f) ? 1.0f : 0.0f;
    const real Htheta_v = (Vm - theta_v > 0.0f) ? 1.0f : 0.0f;

    const real tau_o = (1.0f - Htheta_o) * tau_o1 + Htheta_o * tau_o2;
    const real tau_so = tau_so1 + (tau_so2 - tau_so1) * (1.0f + tanh(k_so * (Vm - u_so))) * 0.5f;

    // Currents
    const real J_fi = -v * Htheta_v * (Vm - theta_v) * (u_u - Vm) / tau_fi;
    const real J_so = ((Vm - u_o) * (1.0f - Htheta_w) / tau_o) + (Htheta_w / tau_so);
    const real J_si = -Htheta_w * w * s / tau_si;

    // Compute the derivative of Vm
    return J_fi + J_so + J_si;
}

static inline void update_sVtilde_MV(real *restrict sVtilde, const real Vm, const real *restrict rhs_sV, const real delta_t)
{
    // Extract state variables
    const real v = rhs_sV[0];
    const real w = rhs_sV[1];
    const real s = rhs_sV[2];

    // Auxiliary variables
    const real Htheta_w = (Vm - theta_w > 0.0f) ? 1.0f : 0.0f;
    const real Htheta_o = (Vm - theta_o > 0.0f) ? 1.0f : 0.0f;
    const real Htheta_v = (Vm - theta_v > 0.0f) ? 1.0f : 0.0f;
    const real Htheta_vminus = (Vm - theta_vminus > 0.0f) ? 1.0f : 0.0f;
    const real tau_vminus = (1.0f - Htheta_vminus) * tau_v1minus + Htheta_vminus * tau_v2minus;

    // Auxiliary variables for Rush-Larsen or Forward Euler
    const real v_inf = ((Vm < theta_vminus) ? 1.0f : 0.0f);
    const real tau_v_RL = (tau_vplus * tau_vminus) / (tau_vplus - tau_vplus * Htheta_v + tau_vminus * Htheta_v);
    const real v_inf_RL = (tau_vplus * v_inf * (1.0f - Htheta_v)) / (tau_vplus - tau_vplus * Htheta_v + tau_vminus * Htheta_v);

    const real w_inf = ((1.0f - Htheta_o) * (1.0f - (Vm / tau_winf)) + Htheta_o * w_infstar);
    const real tau_wminus = tau_w1minus + (tau_w2minus - tau_w1minus) * (1.0f + tanh(k_wminus * (Vm - u_wminus))) * 0.5f;
    const real tau_w_RL = (tau_wplus * tau_wminus) / (tau_wplus - tau_wplus * Htheta_w + tau_wminus * Htheta_w);
    const real w_inf_RL = (tau_wplus * w_inf * (1.0f - Htheta_w)) / (tau_wplus - tau_wplus * Htheta_w + tau_wminus * Htheta_w);

    const real tau_s = (1.0f - Htheta_w) * tau_s1 + Htheta_w * tau_s2;
    const real s_inf_RL = (1.0f + tanh(k_s * (Vm - u_s))) * 0.5f;

    // Calculate approximations with Rush-Larsen or Forward Euler using half time step
    (tau_v_RL > 1.0e-10)
        ? (sVtilde[0] = v_inf_RL - (v_inf_RL - v) * exp(-delta_t / tau_v_RL))
        : (sVtilde[0] = v + delta_t * (1.0f - Htheta_v) * (v_inf - v) / tau_vminus - Htheta_v * v / tau_vplus);

    (tau_w_RL > 1.0e-10)
        ? (sVtilde[1] = w_inf_RL - (w_inf_RL - w) * exp(-delta_t / tau_w_RL))
        : (sVtilde[1] = w + delta_t * (1.0f - Htheta_w) * (w_inf - w) / tau_wminus - Htheta_w * w / tau_wplus);

    (tau_s > 1.0e-10)
        ? (sVtilde[2] = s_inf_RL - (s_inf_RL - s) * exp(-delta_t / tau_s))
        : (sVtilde[2] = s + delta_t * (s_inf_RL - s) / tau_s);
}

static inline void update_sV_MV(real *restrict sV, const real *rhs_sV, const real dSdt_Vm, const real *dSdt_sV, const real delta_t, const int idx)
{
    // Extract state variables 
    const real dSdt_v = dSdt_sV[0];
    const real dSdt_w = dSdt_sV[1];
    const real dSdt_s = dSdt_sV[2]; 

    // Auxiliary variables
    const real Htheta_w = (dSdt_Vm - theta_w > 0.0f) ? 1.0f : 0.0f;
    const real Htheta_o = (dSdt_Vm - theta_o > 0.0f) ? 1.0f : 0.0f;
    const real Htheta_v = (dSdt_Vm - theta_v > 0.0f) ? 1.0f : 0.0f;
    const real Htheta_vminus = (dSdt_Vm - theta_vminus > 0.0f) ? 1.0f : 0.0f;
    const real tau_vminus = (1.0f - Htheta_vminus) * tau_v1minus + Htheta_vminus * tau_v2minus;

    // Auxiliary variables for Rush-Larsen or Forward Euler
    const real v_inf = ((dSdt_Vm < theta_vminus) ? 1.0f : 0.0f);
    const real tau_v_RL = (tau_vplus * tau_vminus) / (tau_vplus - tau_vplus * Htheta_v + tau_vminus * Htheta_v);
    const real v_inf_RL = (tau_vplus * v_inf * (1.0f - Htheta_v)) / (tau_vplus - tau_vplus * Htheta_v + tau_vminus * Htheta_v);

    const real w_inf = ((1.0f - Htheta_o) * (1.0f - (dSdt_Vm / tau_winf)) + Htheta_o * w_infstar);
    const real tau_wminus = tau_w1minus + (tau_w2minus - tau_w1minus) * (1.0f + tanh(k_wminus * (dSdt_Vm - u_wminus))) * 0.5f;
    const real tau_w_RL = (tau_wplus * tau_wminus) / (tau_wplus - tau_wplus * Htheta_w + tau_wminus * Htheta_w);
    const real w_inf_RL = (tau_wplus * w_inf * (1.0f - Htheta_w)) / (tau_wplus - tau_wplus * Htheta_w + tau_wminus * Htheta_w);

    const real tau_s = (1.0f - Htheta_w) * tau_s1 + Htheta_w * tau_s2;
    const real s_inf_RL = (1.0f + tanh(k_s * (dSdt_Vm - u_s))) * 0.5f;

    // Extract state variables
    const real v = rhs_sV[0];
    const real w = rhs_sV[1];
    const real s = rhs_sV[2];

    // Calculate approximations with Rush-Larsen or Forward Euler using half time step
    const int idx_sv = idx * MV_NSV;
    (tau_v_RL > 1.0e-10)
        ? (sV[idx_sv] = v_inf_RL - (v_inf_RL - v) * exp(-delta_t / tau_v_RL))
        : (sV[idx_sv] = v + delta_t * (1.0f - Htheta_v) * (v_inf - dSdt_v) / tau_vminus - Htheta_v * dSdt_v / tau_vplus);

    (tau_w_RL > 1.0e-10)
        ? (sV[idx_sv + 1] = w_inf_RL - (w_inf_RL - w) * exp(-delta_t / tau_w_RL))
        : (sV[idx_sv + 1] = w + delta_t * (1.0f - Htheta_w) * (w_inf - dSdt_w) / tau_wminus - Htheta_w * dSdt_w / tau_wplus);

    (tau_s > 1.0e-10)
        ? (sV[idx_sv + 2] = s_inf_RL - (s_inf_RL - s) * exp(-delta_t / tau_s))
        : (sV[idx_sv + 2] = s + delta_t * (s_inf_RL - dSdt_s) / tau_s);
}

#endif // USE_CUDA

static inline const real get_diffusion_coefficient_MV(const real sigma)
{
    return sigma;
}

// Instantiate the MV model
const CellModelSolver MV_MODEL = {
    .n_state_vars = MV_NSV,
    .activation_thershold = MV_ACTIVATION_THRESHOLD,
    .initialize = initialize_MV,
    .get_actual_sV = get_actual_sV_MV,
    .get_diffusion_coefficient = get_diffusion_coefficient_MV,
    .compute_dVmdt = compute_dVmdt_MV,
    .update_sVtilde = update_sVtilde_MV,
    .update_sV = update_sV_MV,
};