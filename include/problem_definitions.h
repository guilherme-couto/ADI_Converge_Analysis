#ifndef PROBLEM_DEFINITIONS_H
#define PROBLEM_DEFINITIONS_H

#include "core_definitions.h"

// Define problem types
#ifdef LINMONO
#define PROBLEM "LINMONO"
#endif

#ifdef MONODOMAIN
#define PROBLEM "MONODOMAIN"
#endif

#ifdef CABLEEQ
#define PROBLEM "CABLEEQ"
#endif

// Define aux structures for MONODOMAIN
#if defined(MONODOMAIN) || defined(CABLEEQ)

typedef struct {
    real amplitude;
    real begin;
    real duration;
    int xMaxDisc;
    int xMinDisc;
    int yMaxDisc;
    int yMinDisc;
} Stimulus;

#endif // MONODOMAIN || CABLEEQ

#endif // PROBLEM_DEFINITIONS_H