#include "../include/save_functions.h"

// Map of save functions
static const struct
{
    const char *name;
    save_function_t function;
    const char *file_extension;
} save_function_map[] = {
    {"save_as_vtk", save_as_vtk, "vtk"},
    {"save_as_txt", save_as_txt, "txt"},
    {NULL, NULL, NULL}};

// Function to get the appropriate save function based on the name
save_function_t get_save_function(const char *name)
{
    if (name == NULL)
    {
        printf("Error: save function name is NULL\n");
        return NULL;
    }

    for (int i = 0; save_function_map[i].name != NULL; i++)
        if (strcmp(save_function_map[i].name, name) == 0)
            return save_function_map[i].function;
    return NULL;
}

// Function to get the appropriate file extension based on the save function name
const char *get_file_extension(const char *name)
{
    if (name == NULL)
    {
        printf("Error: save function name is NULL\n");
        return NULL;
    }

    for (int i = 0; save_function_map[i].name != NULL; i++)
        if (strcmp(save_function_map[i].name, name) == 0)
            return save_function_map[i].file_extension;
    return NULL;
}

// Implementation of save functions
void save_as_vtk(const char *file_path, real *restrict data, int Nx, int Ny, real delta_x, real delta_y)
{
    FILE *file = fopen(file_path, "w");
    if (file == NULL)
    {
        printf("Error opening file %s\n", file_path);
        exit(1);
    }

    // Write the VTK file - Header for 2D data
    fprintf(file, "# vtk DataFile Version 3.0\n");
    fprintf(file, "Monodomain Data 2D\n");
    fprintf(file, "ASCII\n");
    fprintf(file, "DATASET STRUCTURED_POINTS\n");
    fprintf(file, "DIMENSIONS %d %d 1\n", Nx, Ny);
    fprintf(file, "ORIGIN 0 0 0\n");
    fprintf(file, "SPACING %.6g %.6g 1.0\n", delta_x, delta_y);
    fprintf(file, "POINT_DATA %d\n", Nx * Ny);
    fprintf(file, "SCALARS Vm %s 1\n", REAL_TYPE);
    fprintf(file, "LOOKUP_TABLE default\n");

    for (int index = 0; index < Nx * Ny; index++)
        fprintf(file, "%e\n", data[index]);

    fclose(file);
}

void save_as_txt(const char *file_path, real *restrict data, int Nx, int Ny, real delta_x, real delta_y)
{
    FILE *file = fopen(file_path, "w");
    if (file == NULL)
    {
        printf("Error opening file %s\n", file_path);
        exit(1);
    }

    for (int index = 0; index < Nx * Ny; index++)
        fprintf(file, "%e\n", data[index]);

    fclose(file);
}