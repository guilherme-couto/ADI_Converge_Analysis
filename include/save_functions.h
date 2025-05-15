#ifndef SAVE_FUNCTIONS_H
#define SAVE_FUNCTIONS_H

#include "core_definitions.h"

// Prototypes for saving functions
void save_as_vtk(const char *file_path, real *data, int Nx, int Ny, real delta_x, real delta_y);
void save_as_txt(const char *file_path, real *data, int Nx, int Ny, real delta_x, real delta_y);

typedef void (*save_function_t)(const char *, real *, int, int, real, real);

// Get the appropriate save function based on the file extension
save_function_t get_save_function(const char *name);

#endif // SAVE_FUNCTIONS_H

// typedef void (*CellModelSolver)(const SimulationConfig *, Measurement *, real *);

// static const CellModelSolver solvers[] = {
//     solveMonodomainAFHN,   // CELL_MODEL_AFHN
//     solveMonodomainTT2,    // CELL_MODEL_TT2
//     solveMonodomainMV      // CELL_MODEL_MV
// };

// if (config->cell_model >= 0 && config->cell_model < sizeof(solvers)/sizeof(solvers[0])) {
//     solvers[config->cell_model](config, &measurement, time_array);
// } else {
//     ERRORMSG("Invalid cell model selected.");
//     exit(EXIT_FAILURE);
// }
