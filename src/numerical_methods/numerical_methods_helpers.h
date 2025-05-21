#ifndef NUMERICAL_METHODS_HELPERS_H
#define NUMERICAL_METHODS_HELPERS_H

#include "../../include/core_definitions.h"
#include "../../include/config_parser.h"
#include "../../include/auxfuncs.h"
#include "../../include/logger.h"

#ifdef __cplusplus
extern "C" {
#endif

#if !defined(__CUDACC__)

static int lim(const int num, const int N)
{
    return num == -1 ? 1 : (num == N ? N - 2 : num);
}

// Function to get the stimulus value at a given time and position
static real get_stimulus_value(const real actualTime, const int i, const int j,
                                            const Stimulus *stimuli, const int numberOfStimuli)
{
    for (int k = 0; k < numberOfStimuli; k++)
    {
        const Stimulus *s = &stimuli[k];

        if (j < s->x_discretized.min || j > s->x_discretized.max)
            continue;

        if (i < s->y_discretized.min || i > s->y_discretized.max)
            continue;

        return s->amplitude;
    }
    return 0.0f;
}

// Function to compute the diffusion term in 2D considering isotropic diffusion
static real compute_diffusion_term(const real *Vm, const int i, const int j, const int Nx, const int Ny,
                                                const real diff_coeff, const real phi_x, const real phi_y)
{
    // Get the neighboring values
    const int idx = i * Nx + j;
    const int idx_left = i * Nx + lim(j - 1, Nx);
    const int idx_right = i * Nx + lim(j + 1, Nx);
    const int idx_top = lim(i + 1, Ny) * Nx + j;
    const int idx_bottom = lim(i - 1, Ny) * Nx + j;

    const real Vm_center = Vm[idx];
    const real Vm_left = Vm[idx_left];
    const real Vm_right = Vm[idx_right];
    const real Vm_top = Vm[idx_top];
    const real Vm_bottom = Vm[idx_bottom];

    // Compute the diffusion term
    return diff_coeff * (phi_x * (Vm_left - 2.0f * Vm_center + Vm_right) +
                         phi_y * (Vm_top - 2.0f * Vm_center + Vm_bottom));
}

#else // CUDA specific code

static __device__ int lim(const int num, const int N)
{
    return num == -1 ? 1 : (num == N ? N - 2 : num);
}

// Function to get the stimulus value at a given time and position
static __device__ real get_stimulus_value(const real actualTime, const int i, const int j,
                                                                const Stimulus *stimuli, const int numberOfStimuli)
{
    for (int k = 0; k < numberOfStimuli; k++)
    {
        const Stimulus *s = &stimuli[k];
    
        if (actualTime < s->start_time || actualTime > s->start_time + s->duration)
            continue;

        if (j < s->x_discretized.min || j > s->x_discretized.max)
            continue;

        if (i < s->y_discretized.min || i > s->y_discretized.max)
            continue;

        return s->amplitude;
    }
    return 0.0f;
}

// Function to compute the diffusion term in 2D considering isotropic diffusion
static __device__ real compute_diffusion_term(const real *Vm, const int i, const int j, const int Nx, const int Ny,
                                                                    const real diff_coeff, const real phi_x, const real phi_y)
{
    // Get the neighboring values
    const int idx = i * Nx + j;
    const int idx_left = i * Nx + lim(j - 1, Nx);
    const int idx_right = i * Nx + lim(j + 1, Nx);
    const int idx_top = lim(i + 1, Ny) * Nx + j;
    const int idx_bottom = lim(i - 1, Ny) * Nx + j;

    const real Vm_center = Vm[idx];
    const real Vm_left = Vm[idx_left];
    const real Vm_right = Vm[idx_right];
    const real Vm_top = Vm[idx_top];
    const real Vm_bottom = Vm[idx_bottom];

    // Compute the diffusion term
    const real x_term = Vm_left - 2.0f * Vm_center + Vm_right;
    const real y_term = Vm_top - 2.0f * Vm_center + Vm_bottom;

    return diff_coeff * (phi_x * x_term + phi_y * y_term);
}

#endif // USE_CUDA

// Function to get the active stimuli
static int update_and_get_num_active_stimuli(const real actualTime, const Stimulus *stimuli,
                                                          const int numberOfStimuli, Stimulus *active_stimuli)
{
    int num_active_stimuli = 0;
    for (int k = 0; k < numberOfStimuli; k++)
    {
        const Stimulus *s = &stimuli[k];
        if (actualTime >= s->start_time && actualTime <= s->start_time + s->duration)
            active_stimuli[num_active_stimuli++] = *s;
    }
    return num_active_stimuli;
}

// Function to handle the saving of frames
static void handle_frame_saving(const char *pathToSaveData, const char *file_extension,
                                             const save_function_t save_function, const int timeStepCounter,
                                             const real *Vm, const int Nx, const int Ny,
                                             const real delta_x, const real delta_y, const real actualTime)
{
    static char file_path[MAX_STRING_SIZE];
    snprintf(file_path, MAX_STRING_SIZE, "%s/frames/Vm_%05d.%s", pathToSaveData, timeStepCounter, file_extension);
    save_function(file_path, Vm, Nx, Ny, delta_x, delta_y);
    SUCCESSMSG("Frame at time %.2f ms saved to %s\n", actualTime, file_path);
}

// Function to handle the velocity measurement
static void handle_velocity_measurement(const real Vm_x0, const real Vm_x1, real *t0, real *t1, const real thereshold,
                                                     bool *aux_stim_velocity_flag, bool *stim_velocity_measured,
                                                     const real actualTime, const real x0, const real x1, real *stim_velocity)
{
    if (!*aux_stim_velocity_flag)
    {
        if (Vm_x0 > thereshold)
        {
            *t0 = actualTime;
            *aux_stim_velocity_flag = true;
        }
    }
    else
    {
        if (Vm_x1 > thereshold)
        {
            *t1 = actualTime;
            *stim_velocity = ((x1 - x0) / (*t1 - *t0)) * 1000.0f; // cm/ms
            *stim_velocity_measured = true;
            INFOMSG("Stim velocity (measured from %.2f to %.2f cm) is %.4g cm/s\n", x0, x1, *stim_velocity);
        }
    }
}

#ifdef __cplusplus
}
#endif

#endif // NUMERICAL_METHODS_HELPERS_H