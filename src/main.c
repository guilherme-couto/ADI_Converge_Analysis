#include "../include/config_parser.h"
#include "../include/core_definitions.h"
#include "../include/auxfuncs.h"
#include "../include/cpu_functions.h"
#include "../include/logger.h"
#include <signal.h>

#ifdef USE_CUDA
#include "../include/gpu_functions.h"
#endif // USE_CUDA

SimulationConfig config;
bool config_loaded = false;

void handle_sigint(int sig) {
    ERRORMSG("\nSimulation interrupted by user (SIGINT)\n");

    if (config_loaded) {
        free_simulation_config(&config);
    }

    close_logger();
    exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
    // Handle SIGINT (Ctrl+C)
    signal(SIGINT, handle_sigint);

    // Verify arguments
    if (argc < 2)
    {
        printf("Usage: %s <configuration_file>", argv[0]);
        exit(EXIT_FAILURE);
    }

    // Then proceed with normal loading
    int config_status = load_simulation_config(argv[1], &config);
    if (config_status != 0)
    {
        printf("Failed to load configuration (error code: %d)\n", config_status);
        free_simulation_config(&config);
        return -1;
    }
    config_loaded = true;

    // Create directories
    if (createDirectories(config.output_dir, config.remove_old_files) != 0)
    {
        printf("Error creating directories\n");
        free_simulation_config(&config);
        return -1;
    }

    // Populate the stimuli array
    if (populateStimuli(&config) != 0)
    {
        printf("Failed to populate stimuli\n");
        free_simulation_config(&config);
        return -1;
    }

    // Initialize the logger
    char log_filename[MAX_STRING_SIZE];
    snprintf(log_filename, sizeof(log_filename), "%s/output.log", config.output_dir);
    init_logger(log_filename);
    atexit(close_logger); // Ensure logger is closed on exit

    // Log the simulation header
    log_simulation_header(&config, argv[1]);

    // Save a copy of the configuration file
    saveCopyOfSimulationConfig(argv[1], config.output_dir);

    // Check the execution mode and the equation type to determine the simulation type
    if (config.exec_mode == EXEC_SERIAL)
    {
        if (config.equation_type == EQUATION_MONODOMAIN)
        {
            runMonodomainSimulationSerial(&config);
        }
        else
        {
            ERRORMSG("Unsupported equation type");
            free_simulation_config(&config);
            return -1;
        }
    }

#ifdef USE_CUDA

    else if (config.exec_mode == EXEC_CUDA)
    {
        if (config.equation_type == EQUATION_MONODOMAIN)
        {
            runMonodomainSimulationCUDA(&config);
        }
        else
        {
            ERRORMSG("Unsupported equation type");
            free_simulation_config(&config);
            return -1;
        }
    }

#endif // USE_CUDA

    // Free the configuration
    free_simulation_config(&config);

    return EXIT_SUCCESS;
}
