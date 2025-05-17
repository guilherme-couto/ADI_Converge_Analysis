#ifndef SAVE_FUNCTIONS_H
#define SAVE_FUNCTIONS_H

#include "core_definitions.h"

// Prototypes for saving functions
void save_as_vtk(const char *file_path, real *restrict data, int Nx, int Ny, real delta_x, real delta_y);
void save_as_txt(const char *file_path, real *restrict data, int Nx, int Ny, real delta_x, real delta_y);

typedef void (*save_function_t)(const char *, real *restrict, int, int, real, real);

// Get the appropriate save function based on the file extension
save_function_t get_save_function(const char *name);

// Get the appropriate file extension based on the save function name
const char *get_file_extension(const char *name);

#endif // SAVE_FUNCTIONS_H
