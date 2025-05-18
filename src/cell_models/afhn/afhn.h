#ifndef AFHN_MODEL_H
#define AFHN_MODEL_H

#include "../../../include/core_definitions.h"
#include "../../../include/config_parser.h"
#include "../../../include/auxfuncs.h"
#include "../../../include/logger.h"
#include "../cell_models.h"

// Constants for the AFHN model
#define AFHN_CHI 1.0e3f      // cm^-1
#define AFHN_Cm 1.0e-3f      // mF * cm^-2

// Function prototypes
real forcingTerm(real x, real y, real t, real W, real Lx, real Ly, real sigma);

// Functions of the AFHN solver struct
#define AFHN_NSV 1

extern const real G;      // omega^-1 * cm^-2
extern const real eta1;   // omega^-1 * cm^-1
extern const real eta2; // dimensionless
extern const real eta3;   // dimensionless
extern const real vth;   // mV
extern const real vp;

extern const CellModelSolver AFHN_MODEL;

#endif // AFHN_MODEL_H