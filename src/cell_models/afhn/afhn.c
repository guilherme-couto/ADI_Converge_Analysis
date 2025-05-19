#include "afhn.h"

#ifndef USE_CUDA

// Model parameters - Based on Gerardo_Giorda 2007
// const real sigma = 1.2e-3f; // omega^-1 * cm^-1
// const real chi = 1.0e3f;    // cm^-1
// const real Cm = 1.0e-3f;    // mF * cm^-2

const real G = 1.5f;      // omega^-1 * cm^-2
const real eta1 = 4.4f;   // omega^-1 * cm^-1
const real eta2 = 0.012f; // dimensionless
const real eta3 = 1.0f;   // dimensionless
const real vth = 13.0f;   // mV
const real vp = 100.0f;   // mV

// AFHN solver functions
static inline void initialize_AFHN(real *restrict Vm, real *restrict sV, const int Nx, const int Ny)
{
    // Initial conditions
    const real Vm_init = 0.0f;
    const real W_init = 0.0f;

    // Initialize Vm and W with the initial condition
    for (int idx = 0, idx_sv = 0; idx < Nx * Ny; idx++, idx_sv = idx)
    {
        Vm[idx] = Vm_init;
        sV[idx_sv] = W_init;
    }
}

static inline void get_actual_sV_AFHN(real *restrict actualsV, const real *restrict sV, const int idx)
{
    const int idx_sv = idx * AFHN_NSV;
    actualsV[0] = sV[idx_sv];
}

static inline const real compute_dVmdt_AFHN(const real Vm, const real *restrict sV)
{
    const real W = sV[0];
    return (G * Vm * (1.0f - (Vm / vth)) * (1.0f - (Vm / vp))) + (eta1 * Vm * W);
}

static inline void update_sVtilde_AFHN(real *restrict sVtilde, const real Vm, const real *restrict rhs_sV, const real delta_t)
{
    const real W = rhs_sV[0];
    const real dWdt = (eta2 * ((Vm / vp) - (eta3 * W)));

    sVtilde[0] = W + delta_t * dWdt;
}

static inline void update_sV_AFHN(real *restrict sV, const real *rhs_sV, const real dSdt_Vm, const real *dSdt_sV, const real delta_t, const int idx)
{
    const real W = rhs_sV[0];
    const real dSdt_W = dSdt_sV[0];
    const real dWdt = (eta2 * ((dSdt_Vm / vp) - (eta3 * dSdt_W)));

    const int idx_sv = idx * AFHN_NSV;
    sV[idx_sv] = W + delta_t * dWdt;
}

#endif // USE_CUDA

static inline const real get_diffusion_coefficient_AFHN(const real sigma)
{
    return sigma / (AFHN_Cm * AFHN_CHI);
}

// Instantiate the AFHN model
const CellModelSolver AFHN_MODEL = {
    .n_state_vars = AFHN_NSV,
    .activation_thershold = AFHN_ACTIVATION_THRESHOLD,
    .initialize = initialize_AFHN,
    .get_actual_sV = get_actual_sV_AFHN,
    .get_diffusion_coefficient = get_diffusion_coefficient_AFHN,
    .compute_dVmdt = compute_dVmdt_AFHN,
    .update_sVtilde = update_sVtilde_AFHN,
    .update_sV = update_sV_AFHN,
};

// real forcingTerm(real x, real y, real t, real W, real Lx, real Ly, real sigma)
// {
//     // Calculate coefficients for the ADI method
//     const real diff_coeff = sigma / (AFHN_Cm * AFHN_CHI);
//     real exactVm = (exp(-t)) * cos(_PI * x / Lx) * cos(_PI * y / Ly);
//     real reaction = dVmdt(exactVm, W);
//     return (exactVm * (-(AFHN_CHI * AFHN_Cm) + 2.0f * (sigma / (AFHN_CHI * AFHN_Cm)) * _PI * _PI / (Lx * Ly))) + (AFHN_CHI * reaction);
// }