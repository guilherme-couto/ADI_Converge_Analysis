#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "core_definitions.h"

#if defined(SERIAL) || defined(OPENMP)
const real _pi = 3.14159265358979323846f;
#endif

#ifdef GPU
const __constant__ real _pi = 3.14159265358979323846f;
#endif

#endif // CONSTANTS_H