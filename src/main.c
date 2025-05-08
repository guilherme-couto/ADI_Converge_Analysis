#include "../include/config_parser.h"
#include "../include/logger.h"
#include "../include/core_definitions.h"
#include <stdlib.h>

int main(int argc, char *argv[]) {

    // Initialize logger
    log_init(stderr, "simulation.log", true);
    
    // Verify arguments
    if (argc < 2) {
        LOG_ERROR("Usage: %s <config_file.ini>", argv[0]);
        LOG_FATAL("No configuration file provided");
        log_close();
        exit(EXIT_FAILURE);
    }

    // Log the complete configuration file first
    log_file_content(argv[1], "Simulation configuration file content");
    
    // Then proceed with normal loading
    SimulationConfig config;
    LOG_INFO("Loading configuration from: %s", argv[1]);
    int config_status = load_simulation_config(argv[1], &config);
    if (config_status != 0) {
        LOG_ERROR("Failed to load configuration (error code: %d)", config_status);
        log_close();
        exit(EXIT_FAILURE);
    }
    
    LOG_SUCCESS("Configuration loaded successfully");
    LOG_DEBUG("Configuration details:");
    LOG_DEBUG("  Execution mode: %s", 
             config.exec_mode == GPU ? "GPU" : 
             config.exec_mode == OPENMP ? "OpenMP" : "Serial");
    LOG_DEBUG("  Cell model: %s", 
             config.cell_model == AFHN ? "AFHN" :
             config.cell_model == TT2 ? "TT2" : "MV");
    LOG_DEBUG("  dt=%.4f, dx=%.6f, dy=%.6f", config.dt, config.dx, config.dy);
    
    // Inicialização da simulação
    LOG_INFO("Initializing simulation components");
    
    // ... [seu código de inicialização] ...
    LOG_LINE();  // Ponto de verificação
    
    // Simulação principal
    LOG_INFO("Starting simulation with %s method", config.method);
    
    // [seu código de simulação] ...
    //     LOG_DEBUG("Simulation progress: 25%% completed");
    //     // ... [mais simulação] ...
    //     LOG_DEBUG("Simulation progress: 50%% completed");
        
    //     // Exemplo de verificação periódica
    //     if (check_something_wrong()) {
    //         LOG_ERROR("Numerical instability detected at iteration %d", iteration);
    //         // ... tratamento de erro ...
    //     }
        
    //     // ... [continuação da simulação] ...
        
    // LOG_SUCCESS("Simulation completed successfully");
    // LOG_FATAL("Simulation failed: %s", e.what());
    // log_close();
    // exit(EXIT_FAILURE);
    
    // Finalização
    LOG_INFO("Saving results");
    // ... [código para salvar resultados] ...
    
    LOG_INFO("Cleaning up resources");
    // ... [limpeza] ...
    
    LOG_SUCCESS("Program completed successfully");
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