#include "afhn.h"

#ifdef USE_CUDA

// Model parameters - Based on Gerardo_Giorda 2007
// const __constant__ real sigma = 1.2e-3f; // omega^-1 * cm^-1
// const __constant__ real chi = 1.0e3f;    // cm^-1
// const __constant__ real Cm = 1.0e-3f;    // mF * cm^-2

const __constant__ real G = 1.5f;      // omega^-1 * cm^-2
const __constant__ real eta1 = 4.4f;   // omega^-1 * cm^-1
const __constant__ real eta2 = 0.012f; // dimensionless
const __constant__ real eta3 = 1.0f;   // dimensionless
const __constant__ real vth = 13.0f;   // mV
const __constant__ real vp = 100.0f;   // mV

#endif // USE_CUDA