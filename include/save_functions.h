#ifndef SAVE_FUNCTIONS_H
#define SAVE_FUNCTIONS_H

#include "core_definitions.h"

#ifdef __cplusplus
extern "C" {
#endif

// Prototypes for saving functions
void save_as_vtk(const char *file_path, const real *data, const int Nx, const int Ny, const real delta_x, const real delta_y);
void save_as_txt(const char *file_path, const real *data, const int Nx, const int Ny, const real delta_x, const real delta_y);

typedef void (*save_function_t)(const char *, const real *, const int, const int, const real, const real);

// Get the appropriate save function based on the file extension
save_function_t get_save_function(const char *name);

// Get the appropriate file extension based on the save function name
const char *get_file_extension(const char *name);

#ifdef __cplusplus
}
#endif

#endif // SAVE_FUNCTIONS_H
