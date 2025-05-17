#include "cell_models.h"

// Map of cell models and their corresponding solvers
static const struct
{
    const CellModel model;
    cell_model_solver_t solver;
} cell_model_solver_map[] = {
    {CELL_MODEL_AFHN, solveMonodomainAFHN},
    {CELL_MODEL_TT2, solveMonodomainTT2},
    {CELL_MODEL_MV, solveMonodomainMV},
    {CELL_MODEL_INVALID, NULL}};

cell_model_solver_t get_cell_model_solver(const CellModel *model)
{
    if (model == NULL)
    {
        printf("Error: cell model is NULL\n");
        return NULL;
    }

    for (int i = 0; cell_model_solver_map[i].solver != NULL; i++)
        if (cell_model_solver_map[i].model == *model)
            return cell_model_solver_map[i].solver;
    return NULL;
}