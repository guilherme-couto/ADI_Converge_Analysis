#ifndef MV_MODEL_H
#define MV_MODEL_H

#include "../../../include/core_definitions.h"
#include "../../../include/config_parser.h"
#include "../../../include/auxfuncs.h"
#include "../../../include/logger.h"

// Options: ENDO, M, EPI, PB, TNNP -> default is ENDO
#if !defined(MCELL) && !defined(EPI) && !defined(ENDO) && !defined(PB) && !defined(TNNP)
#define ENDO
#endif

void solveMonodomainMV(const SimulationConfig *config, Measurement *measurement, const real *time_array);


// Initial conditions
// const real u_init = 0.0f;
// const real v_init = 1.0f;
// const real w_init = 1.0f;
// const real s_init = 0.0f;

#endif // MV_MODEL_H