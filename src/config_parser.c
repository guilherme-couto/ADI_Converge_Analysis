#include "../include/config_parser.h"
#include "../external/inih/ini.h"

// Convert string to lowercase for case-insensitive comparison
static void strtolower(char *str)
{
    for (; *str; ++str)
        *str = tolower(*str);
}

static int config_parser_handler(void *user, const char *section, const char *name, const char *value)
{
    SimulationConfig *config = (SimulationConfig *)user;
    char lower_value[MAX_STRING_SIZE];

    // Make case-insensitive comparison easier
    strncpy(lower_value, value, sizeof(lower_value));
    strtolower(lower_value);

#define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0

    // Handle the "simulation" section
    if (MATCH("simulation", "execution_mode"))
    {
        if (strstr(lower_value, "gpu") || strstr(lower_value, "cuda"))
            config->exec_mode = EXEC_CUDA;
        else if (strstr(lower_value, "openmp"))
            config->exec_mode = EXEC_OPENMP;
        else if (strstr(lower_value, "serial"))
            config->exec_mode = EXEC_SERIAL;
        else
            config->exec_mode = EXEC_INVALID;
    }
    else if (MATCH("simulation", "equation_type"))
    {
        if (strstr(lower_value, "monodomain"))
            config->equation_type = EQUATION_MONODOMAIN;
        else
            config->equation_type = EQUATION_INVALID;
    }
    else if (MATCH("simulation", "cell_model"))
    {
        if (strstr(lower_value, "afhn"))
            config->cell_model = CELL_MODEL_AFHN;
        else if (strstr(lower_value, "tt2"))
            config->cell_model = CELL_MODEL_TT2;
        else if (strstr(lower_value, "mv"))
            config->cell_model = CELL_MODEL_MV;
        else
            config->cell_model = CELL_MODEL_INVALID;
    }
    else if (MATCH("simulation", "numerical_method"))
    {
        if (strstr(lower_value, "osadi"))
            config->method = METHOD_OSADI;
        else if (strstr(lower_value, "ssiadi"))
            config->method = METHOD_SSIADI;
        else if (strstr(lower_value, "fe"))
            config->method = METHOD_FE;
        else
            config->method = METHOD_INVALID;
    }
    else if (MATCH("simulation", "save_function"))
    {
        strncpy(config->save_function_name, value, sizeof(config->save_function_name));
        config->save_function = get_save_function(config->save_function_name);
    }
    else if (MATCH("simulation", "theta"))
        config->theta = atof(value);
    else if (MATCH("simulation", "dt"))
        config->dt = atof(value);
    else if (MATCH("simulation", "dx"))
        config->dx = atof(value);
    else if (MATCH("simulation", "dy"))
        config->dy = atof(value);
    else if (MATCH("simulation", "sigma"))
        config->sigma = atof(value);
    else if (MATCH("simulation", "total_time"))
        config->total_time = atof(value);
    else if (MATCH("simulation", "length_x"))
        config->Lx = atof(value);
    else if (MATCH("simulation", "length_y"))
        config->Ly = atof(value);
    else if (MATCH("simulation", "save_rate"))
        config->frame_save_rate = atoi(value);
    else if (MATCH("simulation", "number_of_threads"))
        config->number_of_threads = atoi(value);
    else if (MATCH("simulation", "output_dir"))
        strncpy(config->output_dir, value, sizeof(config->output_dir));
    else if (MATCH("simulation", "remove_old_files"))
        config->remove_old_files = (strstr(lower_value, "true") != NULL);
    else if (MATCH("simulation", "path_to_restore_state_files"))
        strncpy(config->path_to_restore_state_files, value, sizeof(config->path_to_restore_state_files));
    else if (MATCH("simulation", "init_mode"))
        strncpy(config->init_mode, value, sizeof(config->init_mode));
    else if (MATCH("simulation", "shift_state"))
        config->shift_state = (strstr(lower_value, "true") != NULL);
    else if (MATCH("simulation", "save_frames"))
        config->save_frames = (strstr(lower_value, "true") != NULL);
    else if (MATCH("simulation", "save_last_frame"))
        config->save_last_frame = (strstr(lower_value, "true") != NULL);
    else if (MATCH("simulation", "save_last_state"))
        config->save_last_state = (strstr(lower_value, "true") != NULL);
    else if (MATCH("simulation", "measure_velocity"))
        config->measure_velocity = (strstr(lower_value, "true") != NULL);

    // Handle the "stimuli" section
    else if (strncmp(section, "stim_", 5) == 0)
    {
        int stim_index = atoi(section + 5) - 1;

        // Validate index
        if (stim_index < 0)
        {
            ERRORMSG("Invalid stimulus index in section: %s -> Stimulus index must be a positive integer\n", section);
            return 0;
        }

        // Resize the stimuli array if necessary
        if (stim_index >= config->stimulus_count)
        {
            Stimulus *new_stim = realloc(config->stimuli, (stim_index + 1) * sizeof(Stimulus));
            if (!new_stim)
            {
                ERRORMSG("Failed to allocate memory for stimuli\n");
                return 0;
            }
            config->stimuli = new_stim;

            // Initialize new stimulus
            memset(&config->stimuli[stim_index], 0, sizeof(Stimulus));
            config->stimulus_count = stim_index + 1;
        }

        // Set stimulus parameters
        Stimulus *stim = &config->stimuli[stim_index];
        if (strcmp(name, "amplitude") == 0)
            stim->amplitude = atof(value);
        else if (strcmp(name, "begin_time") == 0)
            stim->begin_time = atof(value);
        else if (strcmp(name, "duration") == 0)
            stim->duration = atof(value);
        else if (strcmp(name, "x_min") == 0)
            stim->x_range.min = atof(value);
        else if (strcmp(name, "x_max") == 0)
            stim->x_range.max = atof(value);
        else if (strcmp(name, "y_min") == 0)
            stim->y_range.min = atof(value);
        else if (strcmp(name, "y_max") == 0)
            stim->y_range.max = atof(value);
        else
            WARNINGMSG("Unknown configuration parameter: %s.%s -> IGNORING...\n", section, name);
    }
    return 1;
}

int load_simulation_config(const char *filename, SimulationConfig *config)
{
    DEBUGMSG("Loading config from: %s\n", filename);

    memset(config, 0, sizeof(SimulationConfig));

    // Initialize default values
    config->exec_mode = EXEC_INVALID;
    config->equation_type = EQUATION_INVALID;
    config->cell_model = CELL_MODEL_INVALID;
    config->method = METHOD_INVALID;
    config->theta = -1.0f;
    config->dt = -1.0f;
    config->dx = -1.0f;
    config->dy = -1.0f;
    config->sigma = -1.0f;
    config->total_time = -1.0;
    config->Lx = -1.0f;
    config->Ly = -1.0f;
    config->frame_save_rate = -1;
    config->number_of_threads = -1;
    config->shift_state = false;
    config->save_frames = true;
    config->save_last_frame = true;
    config->save_last_state = false;
    config->measure_velocity = true;
    config->stimuli = NULL;
    config->stimulus_count = -1;
    config->save_function = NULL;
    config->save_function_name[0] = '\0';
    config->output_dir[0] = '\0';
    config->remove_old_files = true;
    config->path_to_restore_state_files[0] = '\0';
    config->init_mode[0] = '\0';

    int result = ini_parse(filename, config_parser_handler, config);
    if (result < 0)
    {
        ERRORMSG("Failed to open config file: %s\n", filename);
        return -1;
    }
    if (result > 0)
    {
        ERRORMSG("Error parsing config file at line %d\n", result);
        return -2;
    }

    if (validate_simulation_config(config) == false)
    {
        return -3;
    }

    // Number of steps
    int M = round(config->total_time / config->dt); // Number of time steps
    int Nx = round(config->Lx / config->dx) + 1;    // Spatial steps in x
    int Ny = round(config->Ly / config->dy) + 1;    // Spatial steps in y
    config->M = M;
    config->Nx = Nx;
    config->Ny = Ny;

    SUCCESSMSG("Configuration loaded successfully\n");

    return 0;
}

bool validate_simulation_config(const SimulationConfig *config)
{
    bool valid = true;

    if (config->exec_mode == EXEC_INVALID)
    {
        ERRORMSG("Invalid execution mode specified\n");
        valid = false;
    }
    if (config->equation_type == EQUATION_INVALID)
    {
        ERRORMSG("Invalid equation type specified\n");
        valid = false;
    }
    if (config->cell_model == CELL_MODEL_INVALID)
    {
        ERRORMSG("Invalid cell model specified\n");
        valid = false;
    }
    if (config->method == METHOD_INVALID)
    {
        ERRORMSG("Invalid numerical method specified\n");
        valid = false;
    }
    if (config->dt <= 0)
    {
        ERRORMSG("Time step (dt) must be positive\n");
        valid = false;
    }
    if (config->dx <= 0 || config->dy <= 0)
    {
        ERRORMSG("Spatial steps (dx, dy) must be positive\n");
        valid = false;
    }
    if (config->sigma <= 0)
    {
        ERRORMSG("Diffusion coefficient (sigma) must be positive\n");
        valid = false;
    }
    if (config->total_time <= 0)
    {
        ERRORMSG("Total time must be positive\n");
        valid = false;
    }
    if (config->Lx <= 0 || config->Ly <= 0)
    {
        ERRORMSG("Simulation domain dimensions (Lx, Ly) must be positive\n");
        valid = false;
    }
    if (config->save_function == NULL && (config->save_last_frame || config->save_frames))
    {
        ERRORMSG("Save function must be specified for saving frames\n");
        valid = false;
    }
    if (config->number_of_threads <= 0 && config->exec_mode == EXEC_OPENMP)
    {
        ERRORMSG("Number of threads must be positive for OpenMP execution\n");
        valid = false;
    }
    if (config->frame_save_rate <= 0 && config->save_frames)
    {
        ERRORMSG("Frame save rate must be positive for saving frames\n");
        valid = false;
    }
    if (config->output_dir[0] == '\0')
    {
        ERRORMSG("Output directory must be specified\n");
        valid = false;
    }
    for (int i = 0; i < config->stimulus_count; i++)
    {
        if (config->stimuli[i].begin_time < 0)
        {
            ERRORMSG("Stimulus begin time must be non-negative\n");
            valid = false;
        }
        if (config->stimuli[i].duration <= 0)
        {
            ERRORMSG("Stimulus duration must be positive\n");
            valid = false;
        }
        if (config->stimuli[i].x_range.min < 0 || config->stimuli[i].x_range.max < 0 ||
            config->stimuli[i].y_range.min < 0 || config->stimuli[i].y_range.max < 0)
        {
            ERRORMSG("Stimulus spatial range must be non-negative\n");
            valid = false;
        }
    }

    return valid;
}

void free_simulation_config(SimulationConfig *config)
{
    if (config->stimuli)
    {
        free(config->stimuli);
        config->stimuli = NULL;
    }
    config->stimulus_count = 0;
}

const char *executionModeToString(ExecutionMode mode)
{
    switch (mode)
    {
    case EXEC_SERIAL:
        return "SERIAL";
    case EXEC_OPENMP:
        return "OPENMP";
    case EXEC_CUDA:
        return "CUDA";
    case EXEC_INVALID:
        return "INVALID";
    default:
        return "UNKNOWN";
    }
}

const char *equationTypeToString(EquationType type)
{
    switch (type)
    {
    case EQUATION_MONODOMAIN:
        return "MONODOMAIN";
    case EQUATION_INVALID:
        return "INVALID";
    default:
        return "UNKNOWN";
    }
}

const char *cellModelToString(CellModel model)
{
    switch (model)
    {
    case CELL_MODEL_AFHN:
        return "AFHN";
    case CELL_MODEL_TT2:
        return "TT2";
    case CELL_MODEL_MV:
        return "MV";
    case CELL_MODEL_INVALID:
        return "INVALID";
    default:
        return "UNKNOWN";
    }
}

const char *numericalMethodToString(NumericalMethod method)
{
    switch (method)
    {
    case METHOD_OSADI:
        return "OSADI";
    case METHOD_SSIADI:
        return "SSIADI";
    case METHOD_FE:
        return "FE";
    case METHOD_INVALID:
        return "INVALID";
    default:
        return "UNKNOWN";
    }
}