#include "cell_models.h"
#include "headers.h"

__device__ void select_get_actual_sV(const CellModel cell_model, real *actualsV, const real *sV, const int idx)
{
    switch (cell_model)
    {
        case CELL_MODEL_AFHN:
            d_get_actual_sV_AFHN(actualsV, sV, idx);
            break;
        case CELL_MODEL_MV:
            d_get_actual_sV_MV(actualsV, sV, idx);
            break;
        default:
            printf("Error: Unsupported cell model\n");
            break;
    }
}

__device__ real select_compute_dVmdt(const CellModel cell_model, const real Vm, const real *sV)
{
    switch (cell_model)
    {
        case CELL_MODEL_AFHN:
            return d_compute_dVmdt_AFHN(Vm, sV);
        case CELL_MODEL_MV:
            return d_compute_dVmdt_MV(Vm, sV);
        default:
            printf("Error: Unsupported cell model\n");
            return 0.0f;
    }
}

__device__ void select_update_sVtilde(const CellModel cell_model, real *sVtilde, const real Vm, const real *rhs_sV, const real delta_t)
{
    switch (cell_model)
    {
        case CELL_MODEL_AFHN:
            d_update_sVtilde_AFHN(sVtilde, Vm, rhs_sV, delta_t);
            break;
        case CELL_MODEL_MV:
            d_update_sVtilde_MV(sVtilde, Vm, rhs_sV, delta_t);
            break;
        default:
            printf("Error: Unsupported cell model\n");
            break;
    }
}

__device__ void select_update_sV(const CellModel cell_model, real *sV, const real *rhs_sV, const real dSdt_Vm, const real *dSdt_sV, const real delta_t, const int idx)
{
    switch (cell_model)
    {
        case CELL_MODEL_AFHN:
            d_update_sV_AFHN(sV, rhs_sV, dSdt_Vm, dSdt_sV, delta_t, idx);
            break;
        case CELL_MODEL_MV:
            d_update_sV_MV(sV, rhs_sV, dSdt_Vm, dSdt_sV, delta_t, idx);
            break;
        default:
            printf("Error: Unsupported cell model\n");
            break;
    }
}
