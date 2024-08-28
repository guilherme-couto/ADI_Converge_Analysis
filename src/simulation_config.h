#ifndef SIMULATION_CONFIG_H
#define SIMULATION_CONFIG_H

#include "include.h"

#ifdef SERIAL
const int L = 1;        // Length of each side (cm)
const real T = 0.1;     // Total time (ms)

#ifdef MONOAFHN
// Stimulation parameters
const real stimStrength = 100.0;

// S1
const real stim1Begin = 0.0;              // S1 start time -> ms
const real stim1Duration = 2.0;           // Stimulation duration -> ms
const real stim1xLimit = 0.2;             // Stimulation x limit -> cm
const real stim1yLimit = 2.0;             // Stimulation y limit -> cm ( = L)

// S2
const real stim2Begin = 120.0;            // S2 start time -> ms
const real stim2Duration = 2.0;           // Stimulation duration -> ms
const real stim2xMax = 1.0;               // Stimulation x max -> cm
const real stim2yMax = 1.0;               // Stimulation y max -> cm
const real stim2xMin = 0.0;               // Stimulation x min -> cm
const real stim2yMin = 0.0;               // Stimulation y min -> cm
#endif // MONOAFHN
#endif // SERIAL

#ifdef GPU
const __constant__ int L = 2;   // Length of each side (cm)
const real T = 1;             // Total time (ms)

#ifdef MONOAFHN
// Save frames parameters
const bool saveFrames = false;
const int numberOfFrames = 100;

// Initial conditions
const real V0 = 0.0;
const real W0 = 0.0;

// Stimulation parameters
const __constant__ int numberOfStimuli = 2;            // Number of stimuli
const real stimuliStrength = 0.0;                    // Stimulation strength -> (amplitude)
const real stimuliDuration = 2.0;                      // Stimulation duration -> ms
const real stimuliBegin[] = {0.0, 120.0};              // Stimuli begin time -> ms
const real stimulixMax[] = {0.2, 1.0};                 // Stimuli x max -> cm
const real stimulixMin[] = {0.0, 0.0};                 // Stimuli x min -> cm
const real stimuliyMax[] = {2.0, 1.0};                 // Stimuli y max -> cm
const real stimuliyMin[] = {0.0, 0.0};                 // Stimuli y min -> cm
#endif // MONOAFHN
#endif // GPU

#endif // SIMULATION_CONFIG_H