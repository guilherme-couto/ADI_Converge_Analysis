#ifndef CELL_MODELS_H
#define CELL_MODELS_H

#include "afhn/afhn.h"
#include "tt2/tt2.h"
#include "minimal_ventricular/mv.h"

typedef void (*cell_model_solver_t)(const SimulationConfig *, Measurement *, const real *);

// Get the appropriate function pointer
cell_model_solver_t get_cell_model_solver(const CellModel *cell_model);

#endif // CELL_MODELS_H