#ifndef CPU_FUNCTIONS_H
#define CPU_FUNCTIONS_H

#include "simulation_config.h"

#ifdef DIFF

real exactSolution(real t, real x, real y)
{
    return (exp(-t))*cos(_pi * x) * cos(_pi * y);
}
real forcingTerm(real x, real y, real t)
{
    return exactSolution(t, x, y) * (-1.0f + 2.0f * (_pi * _pi));
}

#endif // DIFF

#ifdef LINMONO

real exactSolution(real t, real x, real y)
{
    return (exp(-t))*cos(_pi * x / L) * cos(_pi * y / L);
}

real forcingTerm(real x, real y, real t)
{
    return exactSolution(t, x, y) * (-(chi * Cm) + chi * G + 2.0f * (sigma / (chi * Cm)) * _pi * _pi / (L * L));
}

#endif // LINMONO

#if defined(MONODOMAIN)

real exactSolution(real t, real x, real y)
{
    return (exp(-t))*cos(_pi * x / Lx) * cos(_pi * y / Ly);
}

#ifdef AFHN

real RHS_Vm(real Vm, real W)
{
    return (G * Vm * (1.0f - (Vm / vth)) * (1.0f - (Vm / vp))) + (eta1 * Vm * W);
}

real RHS_W(real Vm, real W)
{
    return eta2 * ((Vm / vp) - (eta3 * W));
}

real forcingTerm(real x, real y, real t, real W)
{
    real exactVm = exactSolution(t, x, y);
    real reaction = (G * exactVm * (1.0f - (exactVm / vth)) * (1.0f - (exactVm / vp))) + (eta1 * exactVm * W);
    return (exactVm * (-(chi * Cm) + 2.0f * (sigma / (chi * Cm)) * _pi * _pi / (Lx * Ly))) + (chi * reaction);
}

#endif // AFHN
#endif // MONODOMAIN
#endif // CPU_FUNCTIONS_H