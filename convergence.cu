#ifdef GPU
#include "./src/gpu_methods.h"
#else
#include "./src/cpu_methods.h"
#endif

int main(int argc, char *argv[])
{
    // Parameters
    char *method = NULL;
    real delta_t = 0.0;
    real delta_x = 0.0;
    real delta_y = 0.0;
    real theta = 0.0;

    // Read parameters
    #ifndef CABLEEQ
    if (strcmp(argv[1], "theta-ADI") == 0)
    {
        if (argc != 6)
        {
            ERRORMSG("ERROR: usage is %s {method} {delta_t} {delta_x} {delta_y} {theta}\n", argv[0]);
            return 1;
        }
        method = argv[1];
        delta_t = atof(argv[2]);
        delta_x = atof(argv[3]);
        delta_y = atof(argv[4]);
        theta = atof(argv[5]);
    }
    else
    {
        if (argc != 5)
        {
            ERRORMSG("ERROR: usage is %s {method} {delta_t} {delta_x} {delta_y}\n", argv[0]);
            return 1;
        }
        method = argv[1];
        delta_t = atof(argv[2]);
        delta_x = atof(argv[3]);
        delta_y = atof(argv[4]);
    }
    #else
    if (strcmp(argv[1], "theta-ADI") == 0)
    {
        if (argc != 5)
        {
            ERRORMSG("ERROR: usage is %s {method} {delta_t} {delta_x} {theta}\n", argv[0]);
            return 1;
        }
        method = argv[1];
        delta_t = atof(argv[2]);
        delta_x = atof(argv[3]);
        theta = atof(argv[4]);
    }
    else
    {
        if (argc != 4)
        {
            ERRORMSG("ERROR: usage is %s {method} {delta_t} {delta_x}\n", argv[0]);
            return 1;
        }
        method = argv[1];
        delta_t = atof(argv[2]);
        delta_x = atof(argv[3]);
    }
    #endif // CABLEEQ

    // Configuration overview
    printf("CONFIGURATION:\n");
    printf("Real type: %s\n", REAL_TYPE);
    printf("Execution: %s\n", EXECUTION_TYPE);

    #ifdef CONVERGENCE_ANALYSIS
    printf("CONVERGENCE ANALYSIS\n");
    #ifdef MONODOMAIN
    #define AFHN
    printf("By now, Convergence Analysis for MONODOMAIN is only implemented with AFHN cell model\n");
    #endif // MONODOMAIN
    #endif // CONVERGENCE_ANALYSIS

    #ifdef LINMONO
    printf("Problem: LINMONO -> Adapted monodomain with linear reaction (2D)\n");
    #elif defined(DIFFREAC)
    printf("Problem: DIFFREAC -> Diffusion with linear reaction (2D)\n");
    #elif defined(DIFF)
    printf("Problem: DIFF -> Linear diffusion (2D)\n");
    #elif defined(MONODOMAIN)
    printf("Problem: MONODOMAIN -> Monodomain with %s (2D)\n", CELL_MODEL);
    #elif defined(CABLEEQ)
    printf("Problem: CABLEEQ -> Cable equation with %s (1D)\n", CELL_MODEL);
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
    printf("Total time = %.2f ms\n", totalTime);
    printf("\n");

    // Parameters overview
    printf("RUNNING simulation with parameters:\n");
    printf("Method: %s\n", method);
    if (strcmp(method, "theta-ADI") == 0)
        printf("theta: %.2f\n", theta);
    printf("delta_t: %.5g\n", delta_t);
    printf("delta_x: %.5g\n", delta_x);
    #ifndef CABLEEQ
    printf("delta_y: %.5g\n", delta_y);
    #endif // CABLEEQ
    printf("\n");

    // Call function
    RUNSIMULATION(method, delta_t, delta_x, delta_y, theta);

    printf("\nSIMULATION FINISHED!\n");

    return 0;
}