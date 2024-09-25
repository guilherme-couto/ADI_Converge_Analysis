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

    printf("CONFIGURATION:\n");
    #ifdef USE_DOUBLE
    printf("Real type: DOUBLE\n");
    #else
    printf("Real type: FLOAT\n");
    #endif

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

    // Call function
    #ifdef SERIAL
    printf("Execution: SERIAL\n\n");
    printf("RUNNING simulation with parameters:\n");
    printf("Method: %s\n", method);
    printf("delta_t: %f\n", delta_t);
    printf("delta_x: %f\n", delta_x);
    printf("theta: %f\n", theta);
    printf("\n");
    runSimulation(method, delta_t, delta_x, theta);
    #endif
    #ifdef GPU
    printf("Execution: GPU\n\n");
    printf("RUNNING SIMULATION with parameters:\n");
    printf("Method: %s\n", method);
    printf("delta_t: %f\n", delta_t);
    printf("delta_x: %f\n", delta_x);
    printf("theta: %f\n", theta);
    printf("\n");
    runSimulationGPU(method, delta_t, delta_x, theta);
    #endif

    printf("\nSIMULATION FINISHED!\n");

    return 0;
}