#ifndef CPU_FUNCTIONS_H
#define CPU_FUNCTIONS_H

#include "config_parser.h"
#include "core_definitions.h"
#include "logger.h"
#include "auxfuncs.h"
#include "../src/cell_models/cell_models.h"

void runMonodomainSimulationSerial(const SimulationConfig *config);

#endif // CPU_FUNCTIONS_H