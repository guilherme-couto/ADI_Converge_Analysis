#ifdef GPU
#include "./src/gpu_methods.h"
#else
#include "./src/cpu_methods.h"
#endif

int main(int argc, char *argv[])
{
    // Parameters
    real delta_t = 0.0;
    real delta_x = 0.0;
    real delta_y = 0.0;

    // Read parameters
    #ifndef CABLEEQ
    if (argc != 4)
    {
        ERRORMSG("ERROR: usage is %s {delta_t} {delta_x} {delta_y}\n", argv[0]);
        return 1;
    }
    delta_t = atof(argv[1]);
    delta_x = atof(argv[2]);
    delta_y = atof(argv[3]);
    #else // CABLEEQ
    if (argc != 3)
    {
        ERRORMSG("ERROR: usage is %s {delta_t} {delta_x}\n", argv[0]);
        return 1;
    }
    delta_t = atof(argv[1]);
    delta_x = atof(argv[2]);
    #endif // not CABLEEQ

    // Configuration overview
    printf("CONFIGURATION:\n");
    printf("Real type: %s\n", REAL_TYPE);
    printf("Execution: %s\n", EXECUTION_TYPE);

    #ifdef CONVERGENCE_ANALYSIS_FORCING_TERM
    printf("CONVERGENCE ANALYSIS WITH FORCING TERM\n");
    #ifdef MONODOMAIN
    #define AFHN
    printf("By now, Convergence Analysis with forcing term for MONODOMAIN is only implemented with AFHN cell model\n");
    #endif // MONODOMAIN
    #endif // CONVERGENCE_ANALYSIS_FORCING_TERM

    #ifdef LINMONO
    printf("Problem: LINMONO -> Adapted Monodomain with linear reaction (2D)\n");
    #elif defined(MONODOMAIN)
    printf("Problem: MONODOMAIN -> Monodomain with %s (2D)\n", CELL_MODEL);
    #elif defined(CABLEEQ)
    printf("Problem: CABLEEQ -> Cable equation with %s (1D)\n", CELL_MODEL);
    #ifdef GPU
    printf("Cable equation is only implemented for CPU\n");
    return 1;
    #endif // GPU
    #else
    printf("Problem not yet implemented\n");
    return 1;
    #endif

    #ifdef CABLEEQ
    printf("Cable lenght = %.4g cm\n", Lx);
    #else
    printf("Domain length in x = %.4g cm\n", Lx);
    printf("Domain length in y = %.4g cm\n", Ly);
    #endif // CABLEEQ
    printf("Total time = %.4g ms\n", totalTime);
    printf("\n");

    // Parameters overview
    printf("RUNNING simulation with parameters:\n");
    printf("Method: %s\n", METHOD);
    #ifdef THETA
    printf("theta: %.2f\n", THETA);
    #endif // THETA
    printf("delta_t: %.5g ms\n", delta_t);
    printf("delta_x: %.5g cm (%d um)\n", delta_x, CM_TO_UM(delta_x));
    #ifndef CABLEEQ
    printf("delta_y: %.5g cm (%d um)\n", delta_y, CM_TO_UM(delta_y));
    #endif // CABLEEQ
    printf("\n");

    // Call function
    RUNSIMULATION(delta_t, delta_x, delta_y);

    printf("\n");
    printf("EXECUTION FINISHED!\n\n");

    return 0;
}