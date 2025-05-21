#ifndef MONODOMAIN_SOLVER_H
#define MONODOMAIN_SOLVER_H

#include "config_parser.h"
#include "core_definitions.h"
#include "logger.h"
#include "auxfuncs.h"

#ifdef __cplusplus
extern "C" {
#endif

int runMonodomainSimulationSerial(const SimulationConfig *config);
int runMonodomainSimulationOpenMP(const SimulationConfig *config);
int runMonodomainSimulationCUDA(const SimulationConfig *config);

#ifdef __cplusplus
}
#endif

#endif // MONODOMAIN_SOLVER_H