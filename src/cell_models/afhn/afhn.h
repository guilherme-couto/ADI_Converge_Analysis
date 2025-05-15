#ifndef AFHN_MODEL_H
#define AFHN_MODEL_H

#include "../../../include/core_definitions.h"
#include "../../../include/config_parser.h"
#include "../../../include/auxfuncs.h"

// Function prototypes
real forcingTerm(real x, real y, real t, real W, real Lx, real Ly, real sigma);
void initializeWithInitialConditionAFHN(const real Nx, const real Ny, real *Vm, real *W);

void solveMonodomainAFHN(const SimulationConfig *config, Measurement *measurement, const real *time_array);

void run_SSIADI_AFHN(const SimulationConfig *config, Measurement *measurement, const real *time_array, real *Vm, real *W, real *partRHS);
void run_OSADI_AFHN(const SimulationConfig *config, Measurement *measurement, const real *time_array, real *Vm, real *W, real *partRHS);
void run_FE_AFHN(const SimulationConfig *config, Measurement *measurement, const real *time_array, real *Vm, real *W, real *partRHS);

#endif // AFHN_MODEL_H