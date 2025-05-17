#ifndef NUMERICAL_METHODS_H
#define NUMERICAL_METHODS_H

#include "../../include/core_definitions.h"
#include "../../include/config_parser.h"
#include "numerical_methods_helpers.h"
#include "../cell_models/headers.h"

typedef void (*numerical_method_t)(const SimulationConfig *, Measurement *, const real *restrict, const CellModelSolver *, real *restrict, real *restrict, real *restrict);

void runSSIADI(const SimulationConfig *config, Measurement *measurement, const real *restrict time_array, const CellModelSolver *cell_model_solver, real *restrict Vm, real *restrict sV, real *restrict partRHS);
void runOSADI(const SimulationConfig *config, Measurement *measurement, const real *time_array, const CellModelSolver *cell_model_solver, real *Vm, real *sV, real *partRHS);
void runFE(const SimulationConfig *config, Measurement *measurement, const real *time_array, const CellModelSolver *cell_model_solver, real *Vm, real *sV, real *partRHS);

// Function to get the appropriate numerical method
numerical_method_t get_numerical_method(const NumericalMethod *method);

#endif // NUMERICAL_METHODS_H