#include "../include/config_parser.h"
#include "../include/logger.h"
#include "../include/core_definitions.h"
#include "../include/auxfuncs.h"
#include "../include/cpu_functions.h"

#ifdef USE_CUDA
#include "../include/gpu_functions.h"
#endif // USE_CUDA

int main(int argc, char *argv[])
{
    // Verify arguments
    if (argc < 2)
    {
        INFOMSG("Usage: %s <configuration_file>", argv[0]);
        ERRORMSG("No configuration file provided\n");
        exit(EXIT_FAILURE);
    }

    // Then proceed with normal loading
    SimulationConfig config;
    INFOMSG("Loading configuration from: %s\n", argv[1]);
    int config_status = load_simulation_config(argv[1], &config);
    if (config_status != 0)
    {
        ERRORMSG("Failed to load configuration (error code: %d)\n", config_status);
        exit(EXIT_FAILURE);
    }

    // Create directories
    if (createDirectories(config.output_dir, config.remove_old_files) != 0)
    {
        printf("Error creating directories\n");
        exit(EXIT_FAILURE);
    }
    INFOMSG("Output directory: %s\n", config.output_dir);

    // Populate the stimuli array
    if (populateStimuli(&config) != 0)
    {
        ERRORMSG("Failed to populate stimuli\n");
        free_simulation_config(&config);
        return -4;
    }

    // Print the simulation configuration
    printSimulationConfig(&config);

    // Initialize logger
    char log_file[MAX_STRING_SIZE];
    snprintf(log_file, MAX_STRING_SIZE, "%s/simulation.log", config.output_dir);
    INFOMSG("Log file: %s\n", log_file);

    if (!log_init(stderr, log_file, true))
    {
        ERRORMSG("Failed to initialize logger\n");
        free_simulation_config(&config);
        return -3;
    }
    
    // Log the complete configuration file first
    // log_file_content(argv[1], "Simulation configuration file content");

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
            LOG_ERROR("Unsupported equation type");
            log_close();
            return -2;
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
            LOG_ERROR("Unsupported equation type");
            log_close();
            return -2;
        }
    }

#endif // USE_CUDA


    log_close();
    return EXIT_SUCCESS;
}

// log_init(stderr, "debug.log", true);  // Cores ativas e log detalhado
// log_init(NULL, "production.log", false);  // Sem terminal, só arquivo
// #ifdef DEBUG_MODE
// LOG_DEBUG("State vector: %s", format_vector(state));
// #endif
// // Evite em loops apertados:
// for (int i = 0; i < 1000000; i++) {
//     // NÃO FAÇA: LOG_DEBUG("Iteration %d", i);
//     if (i % 10000 == 0) LOG_DEBUG("Progress: %d iterations", i);
// }