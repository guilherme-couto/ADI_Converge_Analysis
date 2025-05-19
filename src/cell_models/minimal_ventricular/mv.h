#ifndef MV_MODEL_H
#define MV_MODEL_H

#include "../../../include/core_definitions.h"
#include "../cell_models.h"

// Options: ENDO, M, EPI, PB, TNNP -> default is ENDO
#if !defined(MCELL) && !defined(EPI) && !defined(ENDO) && !defined(PB) && !defined(TNNP)
#define ENDO
#endif

#define MV_NSV 3 // Number of state variables in the MV model
#define MV_ACTIVATION_THRESHOLD 0.8f // Activation threshold for the MV model

// Parameters - Based on Minimal Ventricular model
// Model definition https://www.sciencedirect.com/science/article/pii/S0022519308001690?via%3Dihub

// const real Dtildef * 1e-3; // cm^2/s
extern const real u_o;
extern const real u_u;
extern const real theta_v;
extern const real theta_w;
extern const real theta_vminus;
extern const real theta_o;
extern const real tau_v1minus;
extern const real tau_v2minus;
extern const real tau_vplus;
extern const real tau_w1minus;
extern const real tau_w2minus;
extern const real k_wminus;
extern const real u_wminus;
extern const real tau_wplus;
extern const real tau_fi;
extern const real tau_o1;
extern const real tau_o2;
extern const real tau_so1;
extern const real tau_so2;
extern const real k_so;
extern const real u_so;
extern const real tau_s1;
extern const real tau_s2;
extern const real k_s;
extern const real u_s;
extern const real tau_si;
extern const real tau_winf;
extern const real w_infstar;

// Rescale Vm -> from Minimal Ventricular paper
static inline const real rescaleVm(real Vm)
{
    return 85.7f * Vm - 84.0f;
}

extern const CellModelSolver MV_MODEL;

#endif // MV_MODEL_H
