#ifndef CONFIG_PARSER_H
#define CONFIG_PARSER_H

#include <stdbool.h>
#include "core_definitions.h"

// Execution modes
typedef enum {
    SERIAL,
    OPENMP,
    GPU,
    INVALID
} ExecutionMode;

// Cell models
typedef enum {
    AFHN,
    TT2,
    MV,
    INVALID
} CellModel;

// Numerical methods
typedef enum {
    ADI,
    OSADI,
    SSIADI,
    THETA_ADI,
    THETA_RK2,
    FE,
    INVALID
} NumericalMethod;

typedef struct {
    // Basic parameters
    ExecutionMode exec_mode;
    CellModel cell_model;
    NumericalMethod method;
    
    // Numeric parameters
    real theta;
    real dt;
    real dx;
    real dy;
    
    // Boolean flags (using stdbool.h)
    bool shift_state;
    bool save_frames;
    bool save_last_frame;
    bool save_last_state;
    bool measure_velocity;
    
    // String parameters
    char init_mode[20];
} SimulationConfig;

int load_simulation_config(const char *filename, SimulationConfig *config);

bool validate_simulation_config(const SimulationConfig *config);

#endif // CONFIG_PARSER_H