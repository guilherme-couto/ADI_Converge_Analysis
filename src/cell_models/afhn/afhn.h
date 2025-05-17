#ifndef AFHN_MODEL_H
#define AFHN_MODEL_H

#include "../../../include/core_definitions.h"
#include "../../../include/config_parser.h"
#include "../../../include/auxfuncs.h"
#include "../../../include/logger.h"
#include "../simulation_helpers.h"

// Constants for the AFHN model
#define AFHN_CHI 1.0e3f      // cm^-1
#define AFHN_Cm 1.0e-3f      // mF * cm^-2

typedef void (*numerical_method_afhn_t)(const SimulationConfig *, Measurement *, const real *, real *restrict, real *restrict, real *restrict);

// Function to get the appropriate numerical method for AFHN model
numerical_method_afhn_t get_numerical_method_afhn(const NumericalMethod *method);

// Function prototypes
real forcingTerm(real x, real y, real t, real W, real Lx, real Ly, real sigma);
void initializeWithInitialConditionAFHN(const real Nx, const real Ny, real *Vm, real *W);

void solveMonodomainAFHN(const SimulationConfig *config, Measurement *measurement, const real *time_array);

#endif // AFHN_MODEL_H