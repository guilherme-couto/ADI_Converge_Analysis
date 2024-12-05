#ifdef GPU
#include "./src/gpu_methods.h"
#else
#include "./src/cpu_methods.h"
#endif

int main(int argc, char *argv[])
{
    // Parameters
    char *method;
    real delta_t;
    real delta_x;
    real theta;

    // Read parameters
    if (argc != 5)
    {
        printf("ERROR: usage is %s {method} {delta_t} {delta_x} {theta}\n", argv[0]);
        return 1;
    }
    method = argv[1];
    delta_t = atof(argv[2]);
    delta_x = atof(argv[3]);
    theta = atof(argv[4]);

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
    printf("\n");

    // Parameters overview
    printf("RUNNING simulation with parameters:\n");
    printf("Method: %s\n", method);
    printf("delta_t: %.5f\n", delta_t);
    printf("delta_x: %.5f\n", delta_x);
    if (strcmp(method, "theta-ADI") == 0)
        printf("theta: %.2f\n", theta);
    printf("\n");

    // Call function
    #ifdef SERIAL
    runSimulation(method, delta_t, delta_x, theta);
    #endif // SERIAL
    #ifdef GPU
    runSimulationGPU(method, delta_t, delta_x, theta);
    #endif // GPU

    printf("\nSIMULATION FINISHED!\n");

    return 0;
}