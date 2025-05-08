#ifndef NUMERICAL_METHODS_H
#define NUMERICAL_METHODS_H

#include "core_definitions.h"

// Define numerical methods
#ifdef ADI // ADI -> Alternating Direction Implicit
#define METHOD "ADI"
#endif

#ifdef OSADI // OS-ADI -> First Order Operator Splitting ADI
#define METHOD "OS-ADI"
#endif

#ifdef SSIADI // SSI-ADI -> Second Order Semi Implicit ADI
#define METHOD "SSI-ADI"
#endif

#ifdef THETASSIADI // theta-ADI -> theta method with ADI
#define METHOD "theta-SSI-ADI"
#endif

#ifdef THETASSIRK2 // theta-RK2 -> theta method with RK2 (only for CABLEEQ)
#define METHOD "theta-SSI-RK2"
#endif

#ifdef FE // FE -> Forward Euler
#define METHOD "FE"
#endif

#endif // NUMERICAL_METHODS_H