#ifndef CPU_FUNCTIONS_H
#define CPU_FUNCTIONS_H

#include "simulation_config.h"

#ifdef DIFF
real exactSolution(real t, real x, real y)
{
    return (exp(-t)) * cos(_pi*x) * cos(_pi*y);
}
real forcingTerm(real x, real y, real t)
{
    return exactSolution(t, x, y) * (-1.0f+2.0f*(_pi*_pi));
}
#endif // DIFF

#ifdef LINMONO
real exactSolution(real t, real x, real y)
{
    return (exp(-t)) * cos(_pi*x/L) * cos(_pi*y/L);
}

real forcingTerm(real x, real y, real t)
{
    return exactSolution(t, x, y) * (-(chi*Cm) + chi*G + 2.0f*(sigma/(chi*Cm))*_pi*_pi/(L*L));
}
#endif // LINMONO

#if defined(MONODOMAIN) || defined(CABLEEQ)
real exactSolution(real t, real x, real y)
{
    return (exp(-t)) * cos(_pi*x/L) * cos(_pi*y/L);
}

#ifdef AFHN
real RHS_V(real V, real W)
{
    return (G*V*(1.0f-(V/vth)) * (1.0f-(V/vp))) + (eta1*V*W);
}

real RHS_W(real V, real W)
{
    return eta2*((V/vp)-(eta3*W));
}

real forcingTerm(real x, real y, real t, real W)
{
    real exactV = exactSolution(t, x, y);
    real reaction = (G*exactV*(1.0f-(exactV/vth)) * (1.0f-(exactV/vp))) + (eta1*exactV*W);
    return (exactV * (-(chi*Cm) + 2.0f*(sigma/(chi*Cm))*_pi*_pi/(L*L))) + (chi*reaction);
}
#endif // AFHN
#endif // MONODOMAIN || CABLEEQ

#endif // CPU_FUNCTIONS_H