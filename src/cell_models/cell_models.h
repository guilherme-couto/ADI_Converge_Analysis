#ifndef CELL_MODELS_H
#define CELL_MODELS_H

#include "../../include/core_definitions.h"
#include "../../include/config_parser.h"
#include "../../include/auxfuncs.h"
#include "../../include/logger.h"

// Types definition for the cell model solvers
typedef void (*initialize_t)(real *restrict, real *restrict, const int, const int);
typedef void (*get_actual_sV_t)(real *restrict, const real *restrict, const int);
typedef const real (*compute_dVmdt_t)(const real, const real *restrict);
typedef void (*update_sVtilde_t)(real *restrict, const real, const real *restrict, const real);
typedef void (*update_sV_t)(real *restrict, const real *, const real, const real *, const real, const int);

// Structure to hold the cell model solver information
typedef struct
{
    int n_state_vars;
    initialize_t initialize;
    get_actual_sV_t get_actual_sV;
    compute_dVmdt_t compute_dVmdt;
    update_sVtilde_t update_sVtilde;
    update_sV_t update_sV;
} CellModelSolver;

// Instantiate the cell model solvers
extern const CellModelSolver AFHN_MODEL;
extern const CellModelSolver TT2_MODEL;
extern const CellModelSolver MV_MODEL;

// Get the cell model solver based on the cell model type
const CellModelSolver *get_solver_struct(const CellModel *cell_model);

#endif // CELL_MODELS_H