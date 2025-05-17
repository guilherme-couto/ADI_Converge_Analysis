#ifndef CELL_MODELS_H
#define CELL_MODELS_H

#include "../../include/core_definitions.h"
#include "../../include/config_parser.h"
#include "../../include/auxfuncs.h"
#include "../../include/logger.h"

// Types definition for the cell model solvers
typedef void (*initialize_t)(real *restrict , real *restrict , const int, const int);
typedef const real (*compute_dVmdt_t)(const real, const real *restrict, const int);
typedef void (*compute_dSdt_t)(const real, const real *restrict, real *restrict, const int);
typedef void (*update_sV_t)(real *, real *, real *restrict, const real, const int);

// Structure to hold the cell model solver information
typedef struct
{
    int n_state_vars;
    initialize_t initialize;
    compute_dVmdt_t compute_dVmdt;
    compute_dSdt_t compute_dSdt;
    update_sV_t update_sV;
} CellModelSolver;

// Instantiate the cell model solvers
extern const CellModelSolver AFHN_MODEL;
extern const CellModelSolver TT2_MODEL;
extern const CellModelSolver MV_MODEL;

// Get the cell model solver based on the cell model type
const CellModelSolver *get_solver_struct(const CellModel *cell_model);

#endif // CELL_MODELS_H