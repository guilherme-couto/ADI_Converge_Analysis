#include "../include/logger.h"
#include <time.h>

FILE *log_file;

void init_logger(const char *filename) {
    log_file = fopen(filename, "w");
    if (!log_file) {
        fprintf(stderr, "Error opening log file '%s'\n", filename);
        return;
    }
}

void close_logger(void) {
    if (log_file) {
        fclose(log_file);
        log_file = NULL;
    }
}

void log_message(const char *color, const char *prefix, const char *fmt, ...) {
    va_list args;

    // Terminal (with colors)
    printf("%s%s%s", color, prefix, RESET);
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);

    // Log file (without colors)
    if (log_file) {
        fprintf(log_file, "%s", prefix);
        va_start(args, fmt);
        vfprintf(log_file, fmt, args);
        va_end(args);
    }
}

static const char* get_current_timestamp()
{
    static char buffer[32];
    time_t now = time(NULL);
    struct tm *t = localtime(&now);
    strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", t);
    return buffer;
}

void log_simulation_header(const SimulationConfig *config, const char *config_filename)
{
    if (!log_file || !config) return;

    fprintf(log_file, "=================================================================\n");
    fprintf(log_file, "                Monodomain Simulation Log\n");
    fprintf(log_file, "-----------------------------------------------------------------\n");
    fprintf(log_file, " Start Time             : %s\n", get_current_timestamp());
    fprintf(log_file, " Configuration File     : %s\n", config_filename);
    fprintf(log_file, " Output Directory       : %s\n", config->output_dir);
    fprintf(log_file, "\n");

    fprintf(log_file, " Execution Mode         : %s\n", executionModeToString(config->exec_mode));
    fprintf(log_file, " Equation Type          : %s\n", equationTypeToString(config->equation_type));
    fprintf(log_file, " Cell Model             : %s\n", cellModelToString(config->cell_model));
    fprintf(log_file, " Numerical Method       : %s\n", numericalMethodToString(config->method));
    fprintf(log_file, "\n");

    fprintf(log_file, " Grid Size (Nx x Ny)    : %d x %d\n", config->Nx, config->Ny);
    fprintf(log_file, " Domain Size (Lx x Ly)  : %.2f x %.2f\n", config->Lx, config->Ly);
    fprintf(log_file, " Spatial Step (dx, dy)  : %.4g, %.4g\n", config->dx, config->dy);
    fprintf(log_file, " Time Step (dt)         : %.6g\n", config->dt);
    fprintf(log_file, " Total Time             : %.2f\n", config->total_time);
    fprintf(log_file, " Number of Steps (M)    : %d\n", config->M);
    fprintf(log_file, " Sigma                  : %.4f\n", config->sigma);
    fprintf(log_file, " Theta                  : %.4f\n", config->theta);
    fprintf(log_file, "\n");

    fprintf(log_file, " Save Frames            : %s\n", config->save_frames ? "Yes" : "No");
    fprintf(log_file, " Save Last Frame        : %s\n", config->save_last_frame ? "Yes" : "No");
    fprintf(log_file, " Save Last State        : %s\n", config->save_last_state ? "Yes" : "No");
    fprintf(log_file, " Shift State            : %s\n", config->shift_state ? "Yes" : "No");
    fprintf(log_file, " Measure Velocity       : %s\n", config->measure_velocity ? "Yes" : "No");
    fprintf(log_file, " Frame Save Rate        : %d\n", config->frame_save_rate);
    fprintf(log_file, " Number of Threads      : %d\n", config->number_of_threads);
    fprintf(log_file, " Restore State Path     : %s\n", config->path_to_restore_state_files);
    fprintf(log_file, " Initial State Mode     : %s\n", config->init_mode);
    fprintf(log_file, "=================================================================\n\n");

    fflush(log_file);
}