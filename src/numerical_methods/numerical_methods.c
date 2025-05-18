#include "numerical_methods.h"

// Map of numerical methods and their corresponding solvers for AFHN model
static const struct
{
    const NumericalMethod method;
    numerical_method_t solver;
} numerical_method_map[] = {
    {METHOD_SSIADI, runSSIADI},
    {METHOD_OSADI, runOSADI},
    {METHOD_FE, runFE},
    {METHOD_INVALID, NULL}};

numerical_method_t get_numerical_method(const NumericalMethod *method)
{
    if (method == NULL)
    {
        printf("Error: numerical method is NULL\n");
        return NULL;
    }

    for (int i = 0; numerical_method_map[i].solver != NULL; i++)
        if (numerical_method_map[i].method == *method)
            return numerical_method_map[i].solver;
    return NULL;
}