#include "../include/config_parser.h"
#include "../include/logger.h"
#include "../external/inih/ini.h"
#include <string.h>
#include <ctype.h>

// Convert string to lowercase for case-insensitive comparison
static void strtolower(char *str) {
    for (; *str; ++str) *str = tolower(*str);
}

static int ini_handler(void* user, const char* section, const char* name, const char* value) {
    SimulationConfig* config = (SimulationConfig*)user;
    char lower_value[32];
    
    // Make case-insensitive comparison easier
    strncpy(lower_value, value, sizeof(lower_value));
    strtolower(lower_value);

    if (strcmp(section, "simulation") == 0) {
        if (strcmp(name, "execution_mode") == 0) {
            if (strstr(lower_value, "gpu")) config->exec_mode = GPU;
            else if (strstr(lower_value, "openmp")) config->exec_mode = OPENMP;
            else config->exec_mode = SERIAL;
        }
        else if (strcmp(name, "cell_model") == 0) {
            if (strstr(lower_value, "afhn")) config->cell_model = AFHN;
            else if (strstr(lower_value, "tt2")) config->cell_model = TT2;
            else if (strstr(lower_value, "mv")) config->cell_model = MV;
            else config->cell_model = INVALID;
        }
        else if (strcmp(name, "method") == 0) {
            if (strstr(lower_value, "adi")) config->method = ADI;
            else if (strstr(lower_value, "osadi")) config->method = OSADI;
            else if (strstr(lower_value, "ssiadi")) config->method = SSIADI;
            else if (strstr(lower_value, "theta-adi")) config->method = THETA_ADI;
            else if (strstr(lower_value, "theta-rk2")) config->method = THETA_RK2;
            else if (strstr(lower_value, "fe")) config->method = FE;
            else config->method = INVALID;
        }
        else if (strcmp(name, "theta") == 0) {
            config->theta = atof(value);
        }
        else if (strcmp(name, "dt") == 0) {
            config->dt = atof(value);
        }
        else if (strcmp(name, "dx") == 0) {
            config->dx = atof(value);
        }
        else if (strcmp(name, "dy") == 0) {
            config->dy = atof(value);
        }
        else if (strcmp(name, "init_mode") == 0) {
            strncpy(config->init_mode, value, sizeof(config->init_mode));
        }
        else if (strcmp(name, "shift_state") == 0) {
            config->shift_state = (strstr(lower_value, "true") != NULL);
        }
        else if (strcmp(name, "save_frames") == 0) {
            config->save_frames = (strstr(lower_value, "true") != NULL);
        }
        else if (strcmp(name, "save_last_frame") == 0) {
            config->save_last_frame = (strstr(lower_value, "true") != NULL);
        }
        else if (strcmp(name, "save_last_state") == 0) {
            config->save_last_state = (strstr(lower_value, "true") != NULL);
        }
        else if (strcmp(name, "measure_velocity") == 0) {
            config->measure_velocity = (strstr(lower_value, "true") != NULL);
        }
        else {
            LOG_WARN("Unknown configuration parameter: %s.%s", section, name);
            return 0;
        }
    }
    return 1;
}

int load_simulation_config(const char *filename, SimulationConfig *config) {
    LOG_DEBUG("Loading config from: %s", filename);

    memset(config, 0, sizeof(SimulationConfig));
    
    int result = ini_parse(filename, ini_handler, config);
    if (result < 0) {
        LOG_ERROR("Failed to open config file: %s", filename);
        return -1;
    }
    if (result > 0) {
        LOG_ERROR("Error parsing config file at line %d", result);
        return -2;
    }
    
    if (!validate_simulation_config(config)) {
        return -3;
    }

    LOG_SUCCESS("Configuration loaded successfully");
    
    return 0;
}

bool validate_simulation_config(const SimulationConfig *config) {
    bool valid = true;
    
    if (config->dt <= 0) {
        LOG_ERROR("Time step (dt) must be positive");
        valid = false;
    }
    if (config->dx <= 0 || config->dy <= 0) {
        LOG_ERROR("Spatial steps (dx, dy) must be positive");
        valid = false;
    }
    if (config->theta < 0 || config->theta > 1) {
        LOG_ERROR("Theta parameter must be between 0 and 1");
        valid = false;
    }
    if (config->cell_model == INVALID) {
        LOG_ERROR("Invalid cell model specified");
        valid = false;
    }
    if (config->method == INVALID) {
        LOG_ERROR("Invalid numerical method specified");
        valid = false;
    }
    
    return valid;
}