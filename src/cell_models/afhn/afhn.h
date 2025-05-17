#ifndef AFHN_MODEL_H
#define AFHN_MODEL_H

#include "../../../include/core_definitions.h"
#include "../../../include/config_parser.h"
#include "../../../include/auxfuncs.h"
#include "../../../include/logger.h"
#include "../simulation_helpers.h"
#include "../cell_models.h"

// Constants for the AFHN model
#define AFHN_CHI 1.0e3f      // cm^-1
#define AFHN_Cm 1.0e-3f      // mF * cm^-2

typedef void (*numerical_method_afhn_t)(const SimulationConfig *, Measurement *, const real *, real *restrict, real *restrict, real *restrict);

// Function to get the appropriate numerical method for AFHN model
numerical_method_afhn_t get_numerical_method_afhn(const NumericalMethod *method);

// Function prototypes
real forcingTerm(real x, real y, real t, real W, real Lx, real Ly, real sigma);
void initializeWithInitialConditionAFHN(const real Nx, const real Ny, real *Vm, real *W);

// Function to solve the monodomain equation with the AFHN model using the specified numerical method
void solveMonodomainAFHN(const SimulationConfig *config, Measurement *measurement, const real *time_array);

// Functions of the AFHN solver struct
#define AFHN_NSV 1

extern const real G;      // omega^-1 * cm^-2
extern const real eta1;   // omega^-1 * cm^-1
extern const real eta2; // dimensionless
extern const real eta3;   // dimensionless
extern const real vth;   // mV
extern const real vp;

// AFHN solver functions
static inline void initialize_AFHN(real *restrict Vm, real *restrict sV, const int Nx, const int Ny)
{
    // Initial conditions
    const real Vm_init = 0.0f;
    const real W_init = 0.0f;

    // Initialize Vm and W with the initial condition
    for (int idx = 0, idx_sv = 0; idx < Nx * Ny; idx++, idx_sv=idx*AFHN_NSV)
    {
        Vm[idx] = Vm_init;
        sV[idx_sv] = W_init;
    }
}

static inline const real compute_dVmdt_AFHN(const real Vm, const real *restrict sV, const int idx)
{
    real W = sV[idx];
    const real result = (G * Vm * (1.0f - (Vm / vth)) * (1.0f - (Vm / vp))) + (eta1 * Vm * W);
    return result;
}

static inline void compute_dSdt_AFHN(const real Vm, const real *restrict sV, real *restrict dSdt, const int idx)
{
    real W = sV[idx];
    dSdt[0] = (eta2 * ((Vm / vp) - (eta3 * W)));
}

static inline void update_sV_AFHN(real *lhs_sV, real *rhs_sV, real *restrict dSdt, const real delta_t, const int idx)
{
    real W = rhs_sV[idx];
    real dWdt = dSdt[0];
    lhs_sV[idx] = W + delta_t * dWdt;
}

extern const CellModelSolver AFHN_MODEL;

#endif // AFHN_MODEL_H