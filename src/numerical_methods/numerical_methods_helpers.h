#ifndef NUMERICAL_METHODS_HELPERS_H
#define NUMERICAL_METHODS_HELPERS_H

#include "../../include/core_definitions.h"
#include "../../include/config_parser.h"
#include "../../include/auxfuncs.h"
#include "../../include/logger.h"

#ifndef USE_CUDA

// Function to get the stimulus value at a given time and position
static inline real get_stimulus_valueNM(real actualTime, int i, int j, const Stimulus *stimuli, int numberOfStimuli)
{
    for (int k = 0; k < numberOfStimuli; ++k)
    {
        const Stimulus *s = &stimuli[k];

        if (actualTime < s->begin_time || actualTime > s->begin_time + s->duration)
            continue;

        if (j < s->x_discretized.min || j > s->x_discretized.max)
            continue;

        if (i < s->y_discretized.min || i > s->y_discretized.max)
            continue;

        return s->amplitude;
    }
    return 0.0f;
}

// Funtion to get the diffusion term
static inline real get_diffusion_termNM(const real *Vm, int i, int j, int Nx, int Ny, real coeff, real phi_x, real phi_y)
{
    int idx = i * Nx + j;
    int idx_left = i * Nx + lim(j - 1, Nx);
    int idx_right = i * Nx + lim(j + 1, Nx);
    int idx_top = lim(i + 1, Ny) * Nx + j;
    int idx_bottom = lim(i - 1, Ny) * Nx + j;

    return coeff * (phi_x * (Vm[idx_left] - 2.0f * Vm[idx] + Vm[idx_right]) +
                    phi_y * (Vm[idx_top] - 2.0f * Vm[idx] + Vm[idx_bottom]));
}

#endif // USE_CUDA

// Function to handle the saving of frames
static inline void handle_frame_savingNM(const char *pathToSaveData, const char *file_extension,
                                       const save_function_t save_function, int timeStepCounter,
                                       int frameSaveRate, bool saveFrames, real *restrict Vm, int Nx, int Ny,
                                       real delta_x, real delta_y, real actualTime)
{
    if (!saveFrames || (timeStepCounter % frameSaveRate != 0))
        return;

    static char file_path[MAX_STRING_SIZE];
    snprintf(file_path, MAX_STRING_SIZE, "%s/frames/Vm_%05d.%s", pathToSaveData, timeStepCounter, file_extension);
    save_function(file_path, Vm, Nx, Ny, delta_x, delta_y);
    SUCCESSMSG("Frame at time %.2f ms saved to %s\n", actualTime, file_path);
}

// Function to handle the velocity measurement
static inline void handle_velocity_measurementNM(real *Vm, int idx_x0, int idx_x1,
                                               real *t0, real *t1,
                                               bool *aux_stim_velocity_flag, bool *stim_velocity_measured,
                                               real actualTime, real x0, real x1, real *stim_velocity)
{
    // Check if velocity has already been measured
    if (*stim_velocity_measured)
        return;

    if (!*aux_stim_velocity_flag)
    {
        if (Vm[idx_x0] > 10.0f)
        {
            *t0 = actualTime;
            *aux_stim_velocity_flag = true;
        }
    }
    else
    {
        if (Vm[idx_x1] > 10.0f)
        {
            *t1 = actualTime;
            *stim_velocity = ((x1 - x0) / (*t1 - *t0)) * 1000.0f; // cm/ms
            *stim_velocity_measured = true;
            INFOMSG("Stim velocity (measured from %.2f to %.2f cm) is %.4g cm/s\n", x0, x1, *stim_velocity);
        }
    }
}

// Free allocated resources
static inline void freeResourcesNM(real *Vm, real *sV, real *partRHS)
{
    if (Vm != NULL)
        free(Vm);
    if (sV != NULL)
        free(sV);
    if (partRHS != NULL)
        free(partRHS);
}

#endif // NUMERICAL_METHODS_HELPERS_H