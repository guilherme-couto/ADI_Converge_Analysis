#include "cell_models.h"

static const struct
{
    const CellModel model;
    const CellModelSolver *solver_struct;
} cell_model_solver_map[] = {
    {CELL_MODEL_AFHN, &AFHN_MODEL},
    {CELL_MODEL_MV, &MV_MODEL},
    {CELL_MODEL_INVALID, NULL}};

const CellModelSolver *get_solver_struct(const CellModel *cell_model)
{
    if (cell_model == NULL)
    {
        printf("Error: cell model is NULL\n");
        return NULL;
    }
    for (int i = 0; cell_model_solver_map[i].solver_struct != NULL; i++)
        if (cell_model_solver_map[i].model == *cell_model)
            return cell_model_solver_map[i].solver_struct;
    return NULL;
}