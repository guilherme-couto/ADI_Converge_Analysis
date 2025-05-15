#ifndef GPU_FUNCTIONS_H
#define GPU_FUNCTIONS_H

#include "config_parser.h"
#include "core_definitions.h"
#include "logger.h"
#include "auxfuncs.h"
#include "../src/cell_models/cell_models.h"

#ifdef USE_CUDA

#ifdef AFHN

// Kernel to compute the approximate solution of the reaction-diffusion system and update the state variables
__global__ void computeApprox(int Nx, int Ny, real delta_t, real phi_x, real phi_y, real diff_coeff, real actualTime, Stimulus *d_stimuli, real *d_Vm, real *d_partRHS, real *d_W)
{
    // Obtain the index of the thread
    int i = blockIdx.y * blockDim.y + threadIdx.y; // Coordinate y
    int j = blockIdx.x * blockDim.x + threadIdx.x; // Coordinate x

    if (i < Ny && j < Nx)
    {
        // Calculate the index in the 1D array
        int index = i * Nx + j;

        // Stimulation
        real stim = 0.0f;

#pragma unroll
        for (int si = 0; si < numberOfStimuli; si++)
        {
            const Stimulus &stimulus = d_stimuli[si];
            if (actualTime >= stimulus.begin && actualTime <= stimulus.begin + stimulus.duration &&
                j >= stimulus.xMinDisc && j <= stimulus.xMaxDisc &&
                i >= stimulus.yMinDisc && i <= stimulus.yMaxDisc)
            {
                stim = stimulus.amplitude;
                break;
            }
        }

        // Get the actual central value
        real actualVm = d_Vm[index];

        // State variable
        real actualW = d_W[index];

        // Calculate useful denominators
        real denom_vth = 1.0f / vth;
        real denom_vp = 1.0f / vp;
        real denom_Cm_chi = 1.0f / (Cm * chi);

#if defined(SSIADI) || defined(THETASSIADI)

        // Boundary conditions
        real im1Vm = (i - 1 == -1) ? d_Vm[index + Nx] : d_Vm[index - Nx];
        real ip1Vm = (i + 1 == Ny) ? d_Vm[index - Nx] : d_Vm[index + Nx];
        real jm1Vm = (j - 1 == -1) ? d_Vm[index + 1] : d_Vm[index - 1];
        real jp1Vm = (j + 1 == Nx) ? d_Vm[index - 1] : d_Vm[index + 1];

        real diff_term = diff_coeff * (phi_x * (jm1Vm - 2.0f * actualVm + jp1Vm) + phi_y * (im1Vm - 2.0f * actualVm + ip1Vm));

        // Calculate Vmtilde -> utilde = u^n + 0.5 * dt * (A*u^n + R(u^n))
        real RHS_Vm_term = ((G * actualVm * (1.0f - (actualVm * denom_vth)) * (1.0f - (actualVm * denom_vp))) + (eta1 * actualVm * actualW)) * denom_Cm_chi;
        real Vmtilde = actualVm + 0.5f * diff_term + (0.5f * delta_t * (stim - RHS_Vm_term));

        // Calculate approximation for state variables
        real RHS_W_term = eta2 * ((actualVm * denom_vp) - (eta3 * actualW));
        real Wtilde = actualW + (0.5f * delta_t * RHS_W_term);

        // Preparing part of the RHS of the following linear systems
        real RHS_Vmtilde_term = ((G * Vmtilde * (1.0f - (Vmtilde * denom_vth)) * (1.0f - (Vmtilde * denom_vp))) + (eta1 * Vmtilde * Wtilde)) * denom_Cm_chi; // RHS with Vm* and W*
        d_partRHS[index] = delta_t * (stim - RHS_Vmtilde_term);

        // Update state variables with RK2 -> Wn+1 = Wn + dt*R(Vm*, W*)
        real RHS_Wtilde_term = eta2 * ((Vmtilde * denom_vp) - (eta3 * Wtilde));
        d_W[index] = actualW + delta_t * RHS_Wtilde_term;

#endif // SSIADI || THETASSIADI

#ifdef OSADI

        // Calculate part of the RHS of the following linear systems with Forward Euler
        real RHS_Vm_term = ((G * actualVm * (1.0f - (actualVm * denom_vth)) * (1.0f - (actualVm * denom_vp))) + (eta1 * actualVm * actualW)) * denom_Cm_chi;
        d_partRHS[index] = delta_t * (stim - RHS_Vm_term);

        // Update state variables
        real RHS_W_term = eta2 * ((actualVm * denom_vp) - (eta3 * actualW));
        d_W[index] = actualW + delta_t * RHS_W_term; // with Forward Euler -> Wn+1 = Wn + dt*R(Vmn, Wn)

#endif // OSADI

#ifdef FE

        // Boundary conditions
        real im1Vm = (i - 1 == -1) ? d_Vm[index + Nx] : d_Vm[index - Nx];
        real ip1Vm = (i + 1 == Ny) ? d_Vm[index - Nx] : d_Vm[index + Nx];
        real jm1Vm = (j - 1 == -1) ? d_Vm[index + 1] : d_Vm[index - 1];
        real jp1Vm = (j + 1 == Nx) ? d_Vm[index - 1] : d_Vm[index + 1];

        // Calculate the diffusion term
        real diff_term = diff_coeff * (phi_x * (jm1Vm - 2.0f * actualVm + jp1Vm) + phi_y * (im1Vm - 2.0f * actualVm + ip1Vm));

        // Calculate part of the RHS of the following linear systems with Forward Euler
        real RHS_Vm_term = ((G * actualVm * (1.0f - (actualVm * denom_vth)) * (1.0f - (actualVm * denom_vp))) + (eta1 * actualVm * actualW)) * denom_Cm_chi;
        d_partRHS[index] = actualVm + diff_term + delta_t * (stim - RHS_Vm_term);

        // Update state variables with Forward Euler
        real RHS_W_term = eta2 * ((actualVm * denom_vp) - (eta3 * actualW));
        d_W[index] = actualW + delta_t * RHS_W_term;

#endif // FE
    }
}

#endif // AFHN

#ifdef TT2

// Kernel to compute the approximate solution of the reaction-diffusion system and update the state variables
__global__ void computeApprox(int Nx, int Ny, real delta_t, real phi_x, real phi_y, real diff_coeff, real actualTime, Stimulus *d_stimuli, real *d_Vm, real *d_partRHS, real *d_X_r1, real *d_X_r2, real *d_X_s, real *d_m, real *d_h, real *d_j, real *d_d, real *d_f, real *d_f2, real *d_fCaSS, real *d_s, real *d_r, real *d_Ca_i, real *d_Ca_SR, real *d_Ca_SS, real *d_R_prime, real *d_Na_i, real *d_K_i)
{
    // Obtain the index of the thread
    int i = blockIdx.y * blockDim.y + threadIdx.y; // Coordinate y
    int j = blockIdx.x * blockDim.x + threadIdx.x; // Coordinate x

    if (i < Ny && j < Nx)
    {
        // Calculate the index in the 1D array
        int index = i * Nx + j;

        real stim = 0.0f;

#pragma unroll
        for (int si = 0; si < numberOfStimuli; si++)
        {
            const Stimulus &stimulus = d_stimuli[si];
            if (actualTime >= stimulus.begin && actualTime <= stimulus.begin + stimulus.duration &&
                j >= stimulus.xMinDisc && j <= stimulus.xMaxDisc &&
                i >= stimulus.yMinDisc && i <= stimulus.yMaxDisc)
            {
                stim = stimulus.amplitude;
                break;
            }
        }

        // Load central
        real actualVm = d_Vm[index];

        // Boundary conditions
        real im1Vm = (i - 1 == -1) ? d_Vm[index + Nx] : d_Vm[index - Nx];
        real ip1Vm = (i + 1 == Ny) ? d_Vm[index - Nx] : d_Vm[index + Nx];
        real jm1Vm = (j - 1 == -1) ? d_Vm[index + 1] : d_Vm[index - 1];
        real jp1Vm = (j + 1 == Nx) ? d_Vm[index - 1] : d_Vm[index + 1];

        real diff_term = diff_coeff * (phi_x * (jm1Vm - 2.0f * actualVm + jp1Vm) + phi_y * (im1Vm - 2.0f * actualVm + ip1Vm));

        // State variables
        real actualX_r1 = d_X_r1[index];
        real actualX_r2 = d_X_r2[index];
        real actualX_s = d_X_s[index];
        real actualm = d_m[index];
        real actualh = d_h[index];
        real actualj = d_j[index];
        real actuald = d_d[index];
        real actualf = d_f[index];
        real actualf2 = d_f2[index];
        real actualfCaSS = d_fCaSS[index];
        real actuals = d_s[index];
        real actualr = d_r[index];
        real actualCa_i = d_Ca_i[index];
        real actualCa_SR = d_Ca_SR[index];
        real actualCa_SS = d_Ca_SS[index];
        real actualR_prime = d_R_prime[index];
        real actualNa_i = d_Na_i[index];
        real actualK_i = d_K_i[index];

        // Auxiliary variables
        real VmENa = actualVm - (RTONF * log(Na_o / actualNa_i));
        real E_K = (RTONF * log(K_o / actualK_i));
        real VmEK = actualVm - E_K;
        real alpha_K1 = 0.1f / (1.0f + exp(0.06f * (actualVm - E_K - 200.0f)));
        real beta_K1 = (3.0f * exp(0.0002f * (actualVm - E_K + 100.0f)) + exp(0.1f * (actualVm - E_K - 10.0f))) / (1.0f + exp(-0.5f * (actualVm - E_K)));
        real E_Ks = (RTONF * log((K_o + p_KNa * Na_o) / (actualK_i + p_KNa * actualNa_i)));
        real E_Ca = 0.5f * RTONF * log(Ca_o / actualCa_i);

        // Currents
        real INa = G_Na * (actualm * actualm * actualm) * actualh * actualj * VmENa;
        real IbNa = G_bNa * VmENa;
        real IK1 = G_K1 * (alpha_K1 / (alpha_K1 + beta_K1)) * VmEK;
        real Ito = G_to * actualr * actuals * VmEK;
        real IKr = G_Kr * sqrt(K_o / 5.4f) * actualX_r1 * actualX_r2 * VmEK;
        real IKs = G_Ks * actualX_s * actualX_s * (actualVm - E_Ks);
        real ICaL; // !!!
        (actualVm < 15.0f - 1.0e-5f)
            ? (ICaL = G_CaL * actuald * actualf * actualf2 * actualfCaSS * 4.0f * (actualVm - 15.0f) * (F * F) * (0.25f * actualCa_SS * exp(2.0f * (actualVm - 15.0f) * FONRT) - Ca_o) / (R * T * (exp(2.0f * (actualVm - 15.0f) * FONRT) - 1.0f)))
            : (ICaL = G_CaL * actuald * actualf * actualf2 * actualfCaSS * 2.0f * F * (0.25f * actualCa_SS - Ca_o));
        real INaK = ((((p_KNa * K_o) / (K_o + K_mK)) * actualNa_i) / (actualNa_i + K_mNa)) / (1.0f + (0.1245f * exp(((-0.1f) * actualVm * FONRT))) + (0.0353f * exp(((-actualVm) * FONRT))));
        real INaCa; // !!!
        INaCa = (k_NaCa * ((exp((gamma_I_NaCa * actualVm * FONRT)) * (actualNa_i * actualNa_i * actualNa_i) * Ca_o) - (exp(((gamma_I_NaCa - 1.0f) * actualVm * FONRT)) * (Na_o * Na_o * Na_o) * actualCa_i * alpha))) / (((K_mNa_i * K_mNa_i * K_mNa_i) + (Na_o * Na_o * Na_o)) * (K_mCa + Ca_o) * (1.0f + (k_sat * exp(((gamma_I_NaCa)*actualVm * FONRT)))));
        real IpCa = (G_pCa * actualCa_i) / (K_pCa + actualCa_i);
        real IpK = (G_pK * VmEK) / (1.0f + exp((25.0f - actualVm) / 5.98f));
        real IbCa = G_bCa * (actualVm - E_Ca);

        // RHS_Vm at actual time
        real RHS_Vm_term = INa + IbNa + IK1 + Ito + IKr + IKs + ICaL + INaK + INaCa + IpCa + IpK + IbCa;

#if defined(SSIADI) || defined(THETASSIADI)

        // Calculate Vmtilde -> utilde = u^n + 0.5 * dt * (A*u^n + R(u^n))
        real Vmtilde = actualVm + 0.5f * diff_term + (0.5f * delta_t * (stim - RHS_Vm_term));

        // Preparing part of the RHS of the following linear systems
        // Calculate approximation for state variables
        // Rush-Larsen method - auxiliary variables
        real X_r1_inf = 1.0f / (1.0f + exp((-26.0f - actualVm) / 7.0f));
        real alpha_X_r1 = 450.0f / (1.0f + exp((-45.0f - actualVm) / 10.0f));
        real beta_X_r1 = 6.0f / (1.0f + exp((30.0f + actualVm) / 11.5f));

        real X_r2_inf = 1.0f / (1.0f + exp((actualVm + 88.0f) / 24.0f));
        real alpha_X_r2 = 3.0f / (1.0f + exp((-60.0f - actualVm) / 20.0f));
        real beta_X_r2 = 1.12f / (1.0f + exp((actualVm - 60.0f) / 20.0f));

        real X_s_inf = 1.0f / (1.0f + exp((-5.0f - actualVm) / 14.0f));
        real alpha_X_s = 1400.0f / sqrt(1.0f + exp((5.0f - actualVm) / 6.0f));
        real beta_X_s = 1.0f / (1.0f + exp((-35.0f + actualVm) / 15.0f));

        real m_inf = 1.0f / ((1.0f + exp((-56.86 - actualVm) / 9.03f)) * (1.0f + exp((-56.86f - actualVm) / 9.03f)));
        real alpha_m = 1.0f / (1.0f + exp((-60.0f - actualVm) / 5.0f));
        real beta_m = 0.1f / (1.0f + exp((actualVm + 35.0f) / 5.0f)) + (0.1f / (1.0f + exp((actualVm - 50.0f) / 200.0f)));

        real h_inf = 1.0f / ((1.0f + exp((actualVm + 71.55f) / 7.43f)) * (1.0f + exp((actualVm + 71.55f) / 7.43f)));
        real alpha_h;
        (actualVm < -40.0f)
            ? (alpha_h = 0.057f * exp(-(actualVm + 80.0f) / 6.8f))
            : (alpha_h = 0.0f);
        real beta_h;
        (actualVm < -40.0f)
            ? (beta_h = 2.7f * exp(0.079f * actualVm) + 3.1f * 1.0e5f * exp(0.3485f * actualVm))
            : (beta_h = 0.77f / (0.13f * (1.0f + exp((actualVm + 10.66f) / -11.1f))));

        real j_inf = 1.0f / ((1.0f + exp((actualVm + 71.55f) / 7.43f)) * (1.0f + exp((actualVm + 71.55f) / 7.43f)));
        real alpha_j;
        (actualVm < -40.0f)
            ? (alpha_j = ((-25428.0f * exp(0.2444f * actualVm) - (6.948e-6f * exp((-0.04391f) * actualVm))) * (actualVm + 37.78f)) / (1.0f + exp(0.311f * (actualVm + 79.23f))))
            : (alpha_j = 0.0f);
        real beta_j;
        (actualVm < -40.0f)
            ? (beta_j = (0.02424f * exp(-0.01052f * actualVm)) / (1.0f + exp(-0.1378f * (actualVm + 40.14f))))
            : (beta_j = (0.6f * exp(0.057f * actualVm)) / (1.0f + exp(-0.1f * (actualVm + 32.0f))));

        real inf = 1.0f / (1.0f + exp((-8.0f - actualVm) / 7.5f));
        real alpha_d = 1.4f / (1.0f + exp((-35.0f - actualVm) / 13.0f)) + 0.25f;
        real beta_d = 1.4f / (1.0f + exp((actualVm + 5.0f) / 5.0f));
        real gamma_d = 1.0f / (1.0f + exp((50.0f - actualVm) / 20.0f));

        real f_inf = 1.0f / (1.0f + exp((actualVm + 20.0f) / 7.0f));
        real alpha_f = 1102.5f * exp(-(actualVm + 27.0f) * (actualVm + 27.0f) / 225.0f);
        real beta_f = 200.0f / (1.0f + exp((13.0f - actualVm) / 10.0f));
        real gamma_f = 180.0f / (1.0f + exp((actualVm + 30.0f) / 10.0f)) + 20.0f;

        real f2_inf = 0.67f / (1.0f + exp((actualVm + 35.0f) / 7.0f)) + 0.33f;
        real alpha_f2; // !!!
        alpha_f2 = 562.0f * exp(-(actualVm + 27.0f) * (actualVm + 27.0f) / 240.0f);
        real beta_f2 = 31.0f / (1.0f + exp((25.0f - actualVm) / 10.0f));
        real gamma_f2; // !!!
        gamma_f2 = 80.0f / (1.0f + exp((30.0f + actualVm) / 10.0f));

        real fCaSS_inf = 0.6f / (1.0f + (actualCa_SS * actualCa_SS * 400.0f)) + 0.4f;
        real tau_fCaSS = 80.0f / (1.0f + (actualCa_SS * actualCa_SS * 400.0f)) + 2.0f;

#if defined(EPI) || defined(MCELL)

        real s_inf = 1.0f / (1.0f + exp((actualVm + 20.0f) / 5.0f));
        real tau_s = 85.0f * exp(-(actualVm + 45.0f) * (actualVm + 45.0f) / 320.0f) + 5.0f / (1.0f + exp((actualVm - 20.0f) / 5.0f)) + 3.0f;

#endif // EPI || MCELL
#ifdef ENDO

        real s_inf = 1.0f / (1.0f + exp((actualVm + 28.0f) / 5.0f));
        real tau_s = 1000.0f * exp(-(actualVm + 67.0f) * (actualVm + 67.0f) / 1000.0f) + 8.0f;

#endif // ENDO

        real r_inf = 1.0f / (1.0f + exp((20.0f - actualVm) / 6.0f));
        real tau_r = 9.5f * exp(-(actualVm + 40.0f) * (actualVm + 40.0f) / 1800.0f) + 0.8f;

        // Rush-Larsen method - update approximations
        real X_r1tilde = X_r1_inf - (X_r1_inf - actualX_r1) * exp(-(0.5f * delta_t) / (alpha_X_r1 * beta_X_r1));
        real X_r2tilde = X_r2_inf - (X_r2_inf - actualX_r2) * exp(-(0.5f * delta_t) / (alpha_X_r2 * beta_X_r2));
        real X_stilde = X_s_inf - (X_s_inf - actualX_s) * exp(-(0.5f * delta_t) / (alpha_X_s * beta_X_s + 80.0f));
        real mtilde = m_inf - (m_inf - actualm) * exp(-(0.5f * delta_t) / (alpha_m * beta_m));
        real htilde = h_inf - (h_inf - actualh) * exp(-(0.5f * delta_t) * (alpha_h + beta_h));
        real jtilde = j_inf - (j_inf - actualj) * exp(-(0.5f * delta_t) * (alpha_j + beta_j));
        real dtilde = inf - (inf - actuald) * exp(-(0.5f * delta_t) / (alpha_d * beta_d + gamma_d));
        real ftilde = f_inf - (f_inf - actualf) * exp(-(0.5f * delta_t) / (alpha_f + beta_f + gamma_f));
        real f2tilde = f2_inf - (f2_inf - actualf2) * exp(-(0.5f * delta_t) / (alpha_f2 + beta_f2 + gamma_f2));
        real fCaSStilde = fCaSS_inf - (fCaSS_inf - actualfCaSS) * exp(-(0.5f * delta_t) / tau_fCaSS);
        real stilde = s_inf - (s_inf - actuals) * exp(-(0.5f * delta_t) / tau_s);
        real rtilde = r_inf - (r_inf - actualr) * exp(-(0.5f * delta_t) / tau_r);

        // Explicit method - auxiliary variables
        real Ileak = V_leak * (actualCa_SR - actualCa_i);
        real Iup = V_maxup / (1.0f + (K_up * K_up) / (actualCa_i * actualCa_i));
        real k_CaSR = max_SR - (max_SR - min_SR) / (1.0f + (EC / actualCa_SR) * (EC / actualCa_SR));
        real k1 = k1_prime / k_CaSR;
        real k2 = k2_prime * k_CaSR;
        real O = k1 * actualCa_SS * actualCa_SS * actualR_prime / (k3 + k1 * actualCa_SS * actualCa_SS);
        real Irel = V_rel * O * (actualCa_SR - actualCa_SS);
        real Ixfer = V_xfer * (actualCa_SS - actualCa_i);
        real Ca_i_bufC; // !!!
        Ca_i_bufC = 1.0f / (1.0f + bufC * K_bufC / ((K_bufC + actualCa_i) * (K_bufC + actualCa_i)));
        real Ca_SR_bufSR; // !!!
        Ca_SR_bufSR = 1.0f / (1.0f + bufSR * K_bufSR / ((K_bufSR + actualCa_SR) * (K_bufSR + actualCa_SR)));
        real Ca_SS_bufSS; // !!!
        Ca_SS_bufSS = 1.0f / (1.0f + bufSS * K_bufSS / ((K_bufSS + actualCa_SS) * (K_bufSS + actualCa_SS)));

        // Explicit method - RHS of the state variables
        real RHS_R_prime_term = (k4 * (1.0f - actualR_prime)) - k2 * actualCa_SS * actualR_prime;
        real RHS_Ca_i_term = Ca_i_bufC * ((((Ileak - Iup) * V_SR / V_C) + Ixfer) - (((IbCa + IpCa) - 2.0f * INaCa) * Cm / (2.0f * V_C * F)));
        real RHS_Ca_SR_term = Ca_SR_bufSR * (Iup - Ileak - Irel);
        real RHS_Ca_SS_term = Ca_SS_bufSS * (((-ICaL * Cm / (2.0 * V_SS * F)) + (Irel * V_SR / V_SS)) - (Ixfer * V_C / V_SS));
        real RHS_Na_i_term = -((INa + IbNa + 3.0f * INaK + 3.0f * INaCa) * Cm / (V_C * F));
        real RHS_K_i_term = -((IK1 + Ito + IKr + IKs - 2.0f * INaK + IpK + stim) * Cm / (V_C * F));

        // Explicit method - update approximations
        real R_primetilde = actualR_prime + (0.5f * delta_t * RHS_R_prime_term);
        real Ca_itilde = actualCa_i + (0.5f * delta_t * RHS_Ca_i_term);
        real Ca_SRtilde = actualCa_SR + (0.5f * delta_t * RHS_Ca_SR_term);
        real Ca_SStilde = actualCa_SS + (0.5f * delta_t * RHS_Ca_SS_term);
        real Na_itilde = actualNa_i + (0.5f * delta_t * RHS_Na_i_term);
        real K_itilde = actualK_i + (0.5f * delta_t * RHS_K_i_term);

        // Auxiliary variables with Vmtilde
        real VmENatilde = Vmtilde - (RTONF * log(Na_o / Na_itilde));
        real E_Ktilde = (RTONF * log(K_o / K_itilde));
        real VmEKtilde = Vmtilde - E_Ktilde;
        real alpha_K1tilde = 0.1f / (1.0f + exp(0.06f * (Vmtilde - E_Ktilde - 200.0f)));
        real beta_K1tilde = (3.0f * exp(0.0002f * (Vmtilde - E_Ktilde + 100.0f)) + exp(0.1f * (Vmtilde - E_Ktilde - 10.0f))) / (1.0f + exp(-0.5f * (Vmtilde - E_Ktilde)));
        real E_Kstilde = (RTONF * log((K_o + p_KNa * Na_o) / (K_itilde + p_KNa * Na_itilde)));
        real E_Catilde = 0.5f * RTONF * log(Ca_o / Ca_itilde);

        // Currents with Vmtilde
        real INatilde = G_Na * (mtilde * mtilde * mtilde) * htilde * jtilde * VmENatilde;
        real IbNatilde = G_bNa * VmENatilde;
        real IK1tilde = G_K1 * (alpha_K1tilde / (alpha_K1tilde + beta_K1tilde)) * VmEKtilde;
        real Itotilde = G_to * rtilde * stilde * VmEKtilde;
        real IKrtilde = G_Kr * sqrt(K_o / 5.4f) * X_r1tilde * X_r2tilde * VmEKtilde;
        real IKstilde = G_Ks * X_stilde * X_stilde * (Vmtilde - E_Kstilde);
        real ICaLtilde; // !!!
        (Vmtilde < 15.0f - 1.0e-5f)
            ? (ICaLtilde = G_CaL * dtilde * ftilde * f2tilde * fCaSStilde * 4.0f * (Vmtilde - 15.0f) * (F * F) * (0.25f * Ca_SStilde * exp(2.0f * (Vmtilde - 15.0f) * FONRT) - Ca_o) / (R * T * (exp(2.0f * (Vmtilde - 15.0f) * FONRT) - 1.0f)))
            : (ICaLtilde = G_CaL * dtilde * ftilde * f2tilde * fCaSStilde * 2.0f * F * (0.25f * Ca_SStilde - Ca_o));
        real INaKtilde = ((((p_KNa * K_o) / (K_o + K_mK)) * Na_itilde) / (Na_itilde + K_mNa)) / (1.0f + (0.1245f * exp(((-0.1f) * Vmtilde * FONRT))) + (0.0353f * exp(((-Vmtilde) * FONRT))));
        real INaCatilde; // !!!
        INaCatilde = (k_NaCa * ((exp((gamma_I_NaCa * Vmtilde * FONRT)) * (Na_itilde * Na_itilde * Na_itilde) * Ca_o) - (exp(((gamma_I_NaCa - 1.0f) * Vmtilde * FONRT)) * (Na_o * Na_o * Na_o) * Ca_itilde * alpha))) / (((K_mNa_i * K_mNa_i * K_mNa_i) + (Na_o * Na_o * Na_o)) * (K_mCa + Ca_o) * (1.0f + (k_sat * exp(((gamma_I_NaCa)*Vmtilde * FONRT)))));
        real IpCatilde = (G_pCa * Ca_itilde) / (K_pCa + Ca_itilde);
        real IpKtilde = (G_pK * VmEKtilde) / (1.0f + exp((25.0f - Vmtilde) / 5.98f));
        real IbCatilde = G_bCa * (Vmtilde - E_Catilde);

        // RHS of the main equation with Vmtilde
        real RHS_Vmtilde_term = INatilde + IbNatilde + IK1tilde + Itotilde + IKrtilde + IKstilde + ICaLtilde + INaKtilde + INaCatilde + IpCatilde + IpKtilde + IbCatilde;
        d_partRHS[index] = delta_t * (stim - RHS_Vmtilde_term);

        // Update state variables
        // Rush-Larsen method - auxiliary variables
        real X_r1_inftilde = 1.0f / (1.0f + exp((-26.0f - Vmtilde) / 7.0f));
        real alpha_X_r1tilde = 450.0f / (1.0f + exp((-45.0f - Vmtilde) / 10.0f));
        real beta_X_r1tilde = 6.0f / (1.0f + exp((30.0f + Vmtilde) / 11.5f));

        real X_r2_inftilde = 1.0f / (1.0f + exp((Vmtilde + 88.0f) / 24.0f));
        real alpha_X_r2tilde = 3.0f / (1.0f + exp((-60.0f - Vmtilde) / 20.0f));
        real beta_X_r2tilde = 1.12f / (1.0f + exp((Vmtilde - 60.0f) / 20.0f));

        real X_s_inftilde = 1.0f / (1.0f + exp((-5.0f - Vmtilde) / 14.0f));
        real alpha_X_stilde = 1400.0f / sqrt(1.0f + exp((5.0f - Vmtilde) / 6.0f));
        real beta_X_stilde = 1.0f / (1.0f + exp((-35.0f + Vmtilde) / 15.0f));

        real m_inftilde = 1.0f / ((1.0f + exp((-56.86 - Vmtilde) / 9.03f)) * (1.0f + exp((-56.86f - Vmtilde) / 9.03f)));
        real alpha_mtilde = 1.0f / (1.0f + exp((-60.0f - Vmtilde) / 5.0f));
        real beta_mtilde = 0.1f / (1.0f + exp((Vmtilde + 35.0f) / 5.0f)) + (0.1f / (1.0f + exp((Vmtilde - 50.0f) / 200.0f)));

        real h_inftilde = 1.0f / ((1.0f + exp((Vmtilde + 71.55f) / 7.43f)) * (1.0f + exp((Vmtilde + 71.55f) / 7.43f)));
        real alpha_htilde;
        (Vmtilde < -40.0f)
            ? (alpha_htilde = 0.057f * exp(-(Vmtilde + 80.0f) / 6.8f))
            : (alpha_htilde = 0.0f);
        real beta_htilde;
        (Vmtilde < -40.0f)
            ? (beta_htilde = 2.7f * exp(0.079f * Vmtilde) + 3.1f * 1.0e5f * exp(0.3485f * Vmtilde))
            : (beta_htilde = 0.77f / (0.13f * (1.0f + exp((Vmtilde + 10.66f) / -11.1f))));

        real j_inftilde = 1.0f / ((1.0f + exp((Vmtilde + 71.55f) / 7.43f)) * (1.0f + exp((Vmtilde + 71.55f) / 7.43f)));
        real alpha_jtilde;
        (Vmtilde < -40.0f)
            ? (alpha_jtilde = ((-25428.0f * exp(0.2444f * Vmtilde) - (6.948e-6f * exp((-0.04391f) * Vmtilde))) * (Vmtilde + 37.78f)) / (1.0f + exp(0.311f * (Vmtilde + 79.23f))))
            : (alpha_jtilde = 0.0f);
        real beta_jtilde;
        (Vmtilde < -40.0f)
            ? (beta_jtilde = (0.02424f * exp(-0.01052f * Vmtilde)) / (1.0f + exp(-0.1378f * (Vmtilde + 40.14f))))
            : (beta_jtilde = (0.6f * exp(0.057f * Vmtilde)) / (1.0f + exp(-0.1f * (Vmtilde + 32.0f))));

        real d_inftilde = 1.0f / (1.0f + exp((-8.0f - Vmtilde) / 7.5f));
        real alpha_dtilde = 1.4f / (1.0f + exp((-35.0f - Vmtilde) / 13.0f)) + 0.25f;
        real beta_dtilde = 1.4f / (1.0f + exp((Vmtilde + 5.0f) / 5.0f));
        real gamma_dtilde = 1.0f / (1.0f + exp((50.0f - Vmtilde) / 20.0f));

        real f_inftilde = 1.0f / (1.0f + exp((Vmtilde + 20.0f) / 7.0f));
        real alpha_ftilde = 1102.5f * exp(-(Vmtilde + 27.0f) * (Vmtilde + 27.0f) / 225.0f);
        real beta_ftilde = 200.0f / (1.0f + exp((13.0f - Vmtilde) / 10.0f));
        real gamma_ftilde = 180.0f / (1.0f + exp((Vmtilde + 30.0f) / 10.0f)) + 20.0f;

        real f2_inftilde = 0.67f / (1.0f + exp((Vmtilde + 35.0f) / 7.0f)) + 0.33f;
        real alpha_f2tilde; // !!!
        alpha_f2tilde = 562.0f * exp(-(Vmtilde + 27.0f) * (Vmtilde + 27.0f) / 240.0f);
        real beta_f2tilde = 31.0f / (1.0f + exp((25.0f - Vmtilde) / 10.0f));
        real gamma_f2tilde; // !!!
        gamma_f2tilde = 80.0f / (1.0f + exp((30.0f + Vmtilde) / 10.0f));

        real fCaSS_inftilde = 0.6f / (1.0f + (Ca_SStilde * Ca_SStilde * 400.0f)) + 0.4f;
        real tau_fCaSStilde = 80.0f / (1.0f + (Ca_SStilde * Ca_SStilde * 400.0f)) + 2.0f;

#if defined(EPI) || defined(MCELL)

        real s_inftilde = 1.0f / (1.0f + exp((Vmtilde + 20.0f) / 5.0f));
        real tau_stilde = 85.0f * exp(-(Vmtilde + 45.0f) * (Vmtilde + 45.0f) / 320.0f) + 5.0f / (1.0f + exp((Vmtilde - 20.0f) / 5.0f)) + 3.0f;

#endif // EPI || MCELL
#ifdef ENDO

        real s_inftilde = 1.0f / (1.0f + exp((Vmtilde + 28.0f) / 5.0f));
        real tau_stilde = 1000.0f * exp(-(Vmtilde + 67.0f) * (Vmtilde + 67.0f) / 1000.0f) + 8.0f;

#endif // ENDO

        real r_inftilde = 1.0f / (1.0f + exp((20.0f - Vmtilde) / 6.0f));
        real tau_rtilde = 9.5f * exp(-(Vmtilde + 40.0f) * (Vmtilde + 40.0f) / 1800.0f) + 0.8f;

        // Rush-Larsen method - update approximations
        d_X_r1[index] = X_r1_inftilde - (X_r1_inftilde - actualX_r1) * exp(-delta_t / (alpha_X_r1tilde * beta_X_r1tilde));
        d_X_r2[index] = X_r2_inftilde - (X_r2_inftilde - actualX_r2) * exp(-delta_t / (alpha_X_r2tilde * beta_X_r2tilde));
        d_X_s[index] = X_s_inftilde - (X_s_inftilde - actualX_s) * exp(-delta_t / (alpha_X_stilde * beta_X_stilde + 80.0f));
        d_m[index] = m_inftilde - (m_inftilde - actualm) * exp(-delta_t / (alpha_mtilde * beta_mtilde));
        d_h[index] = h_inftilde - (h_inftilde - actualh) * exp(-delta_t * (alpha_htilde + beta_htilde));
        d_j[index] = j_inftilde - (j_inftilde - actualj) * exp(-delta_t * (alpha_jtilde + beta_jtilde));
        d_d[index] = d_inftilde - (d_inftilde - actuald) * exp(-delta_t / (alpha_dtilde * beta_dtilde + gamma_dtilde));
        d_f[index] = f_inftilde - (f_inftilde - actualf) * exp(-delta_t / (alpha_ftilde + beta_ftilde + gamma_ftilde));
        d_f2[index] = f2_inftilde - (f2_inftilde - actualf2) * exp(-delta_t / (alpha_f2tilde + beta_f2tilde + gamma_f2tilde));
        d_fCaSS[index] = fCaSS_inftilde - (fCaSS_inftilde - actualfCaSS) * exp(-delta_t / tau_fCaSStilde);
        d_s[index] = s_inftilde - (s_inftilde - actuals) * exp(-delta_t / tau_stilde);
        d_r[index] = r_inftilde - (r_inftilde - actualr) * exp(-delta_t / tau_rtilde);

        // Explicit method - auxiliary variables
        real Ileaktilde = V_leak * (Ca_SRtilde - Ca_itilde);
        real Iuptilde = V_maxup / (1.0f + (K_up * K_up) / (Ca_itilde * Ca_itilde));
        real k_CaSRtilde = max_SR - (max_SR - min_SR) / (1.0f + (EC / Ca_SRtilde) * (EC / Ca_SRtilde));
        real k1tilde = k1_prime / k_CaSRtilde;
        real k2tilde = k2_prime * k_CaSRtilde;
        real Otilde = k1tilde * Ca_SStilde * Ca_SStilde * R_primetilde / (k3 + k1tilde * Ca_SStilde * Ca_SStilde);
        real Ireltilde = V_rel * Otilde * (Ca_SRtilde - Ca_SStilde);
        real Ixfertilde = V_xfer * (Ca_SStilde - Ca_itilde);
        real Ca_i_bufCtilde; // !!!
        Ca_i_bufCtilde = 1.0f / (1.0f + bufC * K_bufC / ((K_bufC + Ca_itilde) * (K_bufC + Ca_itilde)));
        real Ca_SR_bufSRtilde; // !!!
        Ca_SR_bufSRtilde = 1.0f / (1.0f + bufSR * K_bufSR / ((K_bufSR + Ca_SRtilde) * (K_bufSR + Ca_SRtilde)));
        real Ca_SS_bufSStilde; // !!!
        Ca_SS_bufSStilde = 1.0f / (1.0f + bufSS * K_bufSS / ((K_bufSS + Ca_SStilde) * (K_bufSS + Ca_SStilde)));

        // Explicit method - RHS of the state variables
        real RHS_R_primetilde_term = (k4 * (1.0f - R_primetilde)) - k2tilde * Ca_SStilde * R_primetilde;
        real RHS_Ca_itilde_term = Ca_i_bufCtilde * ((((Ileaktilde - Iuptilde) * V_SR / V_C) + Ixfertilde) - (((IbCatilde + IpCatilde) - 2.0f * INaCatilde) * Cm / (2.0f * V_C * F)));
        real RHS_Ca_SRtilde_term = Ca_SR_bufSRtilde * (Iuptilde - Ileaktilde - Ireltilde);
        real RHS_Ca_SStilde_term = Ca_SS_bufSStilde * (((-ICaLtilde * Cm / (2.0 * V_SS * F)) + (Ireltilde * V_SR / V_SS)) - (Ixfertilde * V_C / V_SS));
        real RHS_Na_itilde_term = -((INatilde + IbNatilde + 3.0f * INaKtilde + 3.0f * INaCatilde) * Cm / (V_C * F));
        real RHS_K_itilde_term = -((IK1tilde + Itotilde + IKrtilde + IKstilde - 2.0f * INaKtilde + IpKtilde + stim) * Cm / (V_C * F));

        // Explicit method - update approximations with RK2 (Heun's method)
        d_R_prime[index] = actualR_prime + (delta_t * RHS_R_primetilde_term);
        d_Ca_i[index] = actualCa_i + (delta_t * RHS_Ca_itilde_term);
        d_Ca_SR[index] = actualCa_SR + (delta_t * RHS_Ca_SRtilde_term);
        d_Ca_SS[index] = actualCa_SS + (delta_t * RHS_Ca_SStilde_term);
        d_Na_i[index] = actualNa_i + (delta_t * RHS_Na_itilde_term);
        d_K_i[index] = actualK_i + (delta_t * RHS_K_itilde_term);

#endif // SSIADI || THETASSIADI

#ifdef OSADI

        d_partRHS[index] = delta_t * (stim - RHS_Vm_term);

        // Preparing part of the RHS of the following linear systems
        // Rush-Larsen method - auxiliary variables
        real X_r1_inf = 1.0f / (1.0f + exp((-26.0f - actualVm) / 7.0f));
        real alpha_X_r1 = 450.0f / (1.0f + exp((-45.0f - actualVm) / 10.0f));
        real beta_X_r1 = 6.0f / (1.0f + exp((30.0f + actualVm) / 11.5f));

        real X_r2_inf = 1.0f / (1.0f + exp((actualVm + 88.0f) / 24.0f));
        real alpha_X_r2 = 3.0f / (1.0f + exp((-60.0f - actualVm) / 20.0f));
        real beta_X_r2 = 1.12f / (1.0f + exp((actualVm - 60.0f) / 20.0f));

        real X_s_inf = 1.0f / (1.0f + exp((-5.0f - actualVm) / 14.0f));
        real alpha_X_s = 1400.0f / sqrt(1.0f + exp((5.0f - actualVm) / 6.0f));
        real beta_X_s = 1.0f / (1.0f + exp((-35.0f + actualVm) / 15.0f));

        real m_inf = 1.0f / ((1.0f + exp((-56.86 - actualVm) / 9.03f)) * (1.0f + exp((-56.86f - actualVm) / 9.03f)));
        real alpha_m = 1.0f / (1.0f + exp((-60.0f - actualVm) / 5.0f));
        real beta_m = 0.1f / (1.0f + exp((actualVm + 35.0f) / 5.0f)) + (0.1f / (1.0f + exp((actualVm - 50.0f) / 200.0f)));

        real h_inf = 1.0f / ((1.0f + exp((actualVm + 71.55f) / 7.43f)) * (1.0f + exp((actualVm + 71.55f) / 7.43f)));
        real alpha_h;
        (actualVm < -40.0f)
            ? (alpha_h = 0.057f * exp(-(actualVm + 80.0f) / 6.8f))
            : (alpha_h = 0.0f);
        real beta_h;
        (actualVm < -40.0f)
            ? (beta_h = 2.7f * exp(0.079f * actualVm) + 3.1f * 1.0e5f * exp(0.3485f * actualVm))
            : (beta_h = 0.77f / (0.13f * (1.0f + exp((actualVm + 10.66f) / -11.1f))));

        real j_inf = 1.0f / ((1.0f + exp((actualVm + 71.55f) / 7.43f)) * (1.0f + exp((actualVm + 71.55f) / 7.43f)));
        real alpha_j;
        (actualVm < -40.0f)
            ? (alpha_j = ((-25428.0f * exp(0.2444f * actualVm) - (6.948e-6f * exp((-0.04391f) * actualVm))) * (actualVm + 37.78f)) / (1.0f + exp(0.311f * (actualVm + 79.23f))))
            : (alpha_j = 0.0f);
        real beta_j;
        (actualVm < -40.0f)
            ? (beta_j = (0.02424f * exp(-0.01052f * actualVm)) / (1.0f + exp(-0.1378f * (actualVm + 40.14f))))
            : (beta_j = (0.6f * exp(0.057f * actualVm)) / (1.0f + exp(-0.1f * (actualVm + 32.0f))));

        real inf = 1.0f / (1.0f + exp((-8.0f - actualVm) / 7.5f));
        real alpha_d = 1.4f / (1.0f + exp((-35.0f - actualVm) / 13.0f)) + 0.25f;
        real beta_d = 1.4f / (1.0f + exp((actualVm + 5.0f) / 5.0f));
        real gamma_d = 1.0f / (1.0f + exp((50.0f - actualVm) / 20.0f));

        real f_inf = 1.0f / (1.0f + exp((actualVm + 20.0f) / 7.0f));
        real alpha_f = 1102.5f * exp(-(actualVm + 27.0f) * (actualVm + 27.0f) / 225.0f);
        real beta_f = 200.0f / (1.0f + exp((13.0f - actualVm) / 10.0f));
        real gamma_f = 180.0f / (1.0f + exp((actualVm + 30.0f) / 10.0f)) + 20.0f;

        real f2_inf = 0.67f / (1.0f + exp((actualVm + 35.0f) / 7.0f)) + 0.33f;
        real alpha_f2; // !!!
        alpha_f2 = 562.0f * exp(-(actualVm + 27.0f) * (actualVm + 27.0f) / 240.0f);
        real beta_f2 = 31.0f / (1.0f + exp((25.0f - actualVm) / 10.0f));
        real gamma_f2; // !!!
        gamma_f2 = 80.0f / (1.0f + exp((30.0f + actualVm) / 10.0f));

        real fCaSS_inf = 0.6f / (1.0f + (actualCa_SS * actualCa_SS * 400.0f)) + 0.4f;
        real tau_fCaSS = 80.0f / (1.0f + (actualCa_SS * actualCa_SS * 400.0f)) + 2.0f;

#if defined(EPI) || defined(MCELL)

        real s_inf = 1.0f / (1.0f + exp((actualVm + 20.0f) / 5.0f));
        real tau_s = 85.0f * exp(-(actualVm + 45.0f) * (actualVm + 45.0f) / 320.0f) + 5.0f / (1.0f + exp((actualVm - 20.0f) / 5.0f)) + 3.0f;

#endif // EPI || MCELL
#ifdef ENDO

        real s_inf = 1.0f / (1.0f + exp((actualVm + 28.0f) / 5.0f));
        real tau_s = 1000.0f * exp(-(actualVm + 67.0f) * (actualVm + 67.0f) / 1000.0f) + 8.0f;

#endif // ENDO

        real r_inf = 1.0f / (1.0f + exp((20.0f - actualVm) / 6.0f));
        real tau_r = 9.5f * exp(-(actualVm + 40.0f) * (actualVm + 40.0f) / 1800.0f) + 0.8f;

        // Rush-Larsen method - update approximations
        d_X_r1[index] = X_r1_inf - (X_r1_inf - actualX_r1) * exp(-delta_t / (alpha_X_r1 * beta_X_r1));
        d_X_r2[index] = X_r2_inf - (X_r2_inf - actualX_r2) * exp(-delta_t / (alpha_X_r2 * beta_X_r2));
        d_X_s[index] = X_s_inf - (X_s_inf - actualX_s) * exp(-delta_t / (alpha_X_s * beta_X_s + 80.0f));
        d_m[index] = m_inf - (m_inf - actualm) * exp(-delta_t / (alpha_m * beta_m));
        d_h[index] = h_inf - (h_inf - actualh) * exp(-delta_t * (alpha_h + beta_h));
        d_j[index] = j_inf - (j_inf - actualj) * exp(-delta_t * (alpha_j + beta_j));
        d_d[index] = inf - (inf - actuald) * exp(-delta_t / (alpha_d * beta_d + gamma_d));
        d_f[index] = f_inf - (f_inf - actualf) * exp(-delta_t / (alpha_f + beta_f + gamma_f));
        d_f2[index] = f2_inf - (f2_inf - actualf2) * exp(-delta_t / (alpha_f2 + beta_f2 + gamma_f2));
        d_fCaSS[index] = fCaSS_inf - (fCaSS_inf - actualfCaSS) * exp(-delta_t / tau_fCaSS);
        d_s[index] = s_inf - (s_inf - actuals) * exp(-delta_t / tau_s);
        d_r[index] = r_inf - (r_inf - actualr) * exp(-delta_t / tau_r);

        // Explicit method - auxiliary variables
        real Ileak = V_leak * (actualCa_SR - actualCa_i);
        real Iup = V_maxup / (1.0f + (K_up * K_up) / (actualCa_i * actualCa_i));
        real k_CaSR = max_SR - (max_SR - min_SR) / (1.0f + (EC / actualCa_SR) * (EC / actualCa_SR));
        real k1 = k1_prime / k_CaSR;
        real k2 = k2_prime * k_CaSR;
        real O = k1 * actualCa_SS * actualCa_SS * actualR_prime / (k3 + k1 * actualCa_SS * actualCa_SS);
        real Irel = V_rel * O * (actualCa_SR - actualCa_SS);
        real Ixfer = V_xfer * (actualCa_SS - actualCa_i);
        real Ca_i_bufC; // !!!
        Ca_i_bufC = 1.0f / (1.0f + bufC * K_bufC / ((K_bufC + actualCa_i) * (K_bufC + actualCa_i)));
        real Ca_SR_bufSR; // !!!
        Ca_SR_bufSR = 1.0f / (1.0f + bufSR * K_bufSR / ((K_bufSR + actualCa_SR) * (K_bufSR + actualCa_SR)));
        real Ca_SS_bufSS; // !!!
        Ca_SS_bufSS = 1.0f / (1.0f + bufSS * K_bufSS / ((K_bufSS + actualCa_SS) * (K_bufSS + actualCa_SS)));

        // Explicit method - RHS of the state variables
        real RHS_R_prime_term = (k4 * (1.0f - actualR_prime)) - k2 * actualCa_SS * actualR_prime;
        real RHS_Ca_i_term = Ca_i_bufC * ((((Ileak - Iup) * V_SR / V_C) + Ixfer) - (((IbCa + IpCa) - 2.0f * INaCa) * Cm / (2.0f * V_C * F)));
        real RHS_Ca_SR_term = Ca_SR_bufSR * (Iup - Ileak - Irel);
        real RHS_Ca_SS_term = Ca_SS_bufSS * (((-ICaL * Cm / (2.0 * V_SS * F)) + (Irel * V_SR / V_SS)) - (Ixfer * V_C / V_SS));
        real RHS_Na_i_term = -((INa + IbNa + 3.0f * INaK + 3.0f * INaCa) * Cm / (V_C * F));
        real RHS_K_i_term = -((IK1 + Ito + IKr + IKs - 2.0f * INaK + IpK + stim) * Cm / (V_C * F));

        // Explicit method - update approximations
        d_R_prime[index] = actualR_prime + (delta_t * RHS_R_prime_term);
        d_Ca_i[index] = actualCa_i + (delta_t * RHS_Ca_i_term);
        d_Ca_SR[index] = actualCa_SR + (delta_t * RHS_Ca_SR_term);
        d_Ca_SS[index] = actualCa_SS + (delta_t * RHS_Ca_SS_term);
        d_Na_i[index] = actualNa_i + (delta_t * RHS_Na_i_term);
        d_K_i[index] = actualK_i + (delta_t * RHS_K_i_term);

#endif // OSADI

#ifdef FE

        d_partRHS[index] = actualVm + diff_term + (delta_t * (stim - RHS_Vm_term));

        // Preparing part of the RHS of the following linear systems
        // Rush-Larsen method - auxiliary variables
        real X_r1_inf = 1.0f / (1.0f + exp((-26.0f - actualVm) / 7.0f));
        real alpha_X_r1 = 450.0f / (1.0f + exp((-45.0f - actualVm) / 10.0f));
        real beta_X_r1 = 6.0f / (1.0f + exp((30.0f + actualVm) / 11.5f));

        real X_r2_inf = 1.0f / (1.0f + exp((actualVm + 88.0f) / 24.0f));
        real alpha_X_r2 = 3.0f / (1.0f + exp((-60.0f - actualVm) / 20.0f));
        real beta_X_r2 = 1.12f / (1.0f + exp((actualVm - 60.0f) / 20.0f));

        real X_s_inf = 1.0f / (1.0f + exp((-5.0f - actualVm) / 14.0f));
        real alpha_X_s = 1400.0f / sqrt(1.0f + exp((5.0f - actualVm) / 6.0f));
        real beta_X_s = 1.0f / (1.0f + exp((-35.0f + actualVm) / 15.0f));

        real m_inf = 1.0f / ((1.0f + exp((-56.86 - actualVm) / 9.03f)) * (1.0f + exp((-56.86f - actualVm) / 9.03f)));
        real alpha_m = 1.0f / (1.0f + exp((-60.0f - actualVm) / 5.0f));
        real beta_m = 0.1f / (1.0f + exp((actualVm + 35.0f) / 5.0f)) + (0.1f / (1.0f + exp((actualVm - 50.0f) / 200.0f)));

        real h_inf = 1.0f / ((1.0f + exp((actualVm + 71.55f) / 7.43f)) * (1.0f + exp((actualVm + 71.55f) / 7.43f)));
        real alpha_h;
        (actualVm < -40.0f)
            ? (alpha_h = 0.057f * exp(-(actualVm + 80.0f) / 6.8f))
            : (alpha_h = 0.0f);
        real beta_h;
        (actualVm < -40.0f)
            ? (beta_h = 2.7f * exp(0.079f * actualVm) + 3.1f * 1.0e5f * exp(0.3485f * actualVm))
            : (beta_h = 0.77f / (0.13f * (1.0f + exp((actualVm + 10.66f) / -11.1f))));

        real j_inf = 1.0f / ((1.0f + exp((actualVm + 71.55f) / 7.43f)) * (1.0f + exp((actualVm + 71.55f) / 7.43f)));
        real alpha_j;
        (actualVm < -40.0f)
            ? (alpha_j = ((-25428.0f * exp(0.2444f * actualVm) - (6.948e-6f * exp((-0.04391f) * actualVm))) * (actualVm + 37.78f)) / (1.0f + exp(0.311f * (actualVm + 79.23f))))
            : (alpha_j = 0.0f);
        real beta_j;
        (actualVm < -40.0f)
            ? (beta_j = (0.02424f * exp(-0.01052f * actualVm)) / (1.0f + exp(-0.1378f * (actualVm + 40.14f))))
            : (beta_j = (0.6f * exp(0.057f * actualVm)) / (1.0f + exp(-0.1f * (actualVm + 32.0f))));

        real inf = 1.0f / (1.0f + exp((-8.0f - actualVm) / 7.5f));
        real alpha_d = 1.4f / (1.0f + exp((-35.0f - actualVm) / 13.0f)) + 0.25f;
        real beta_d = 1.4f / (1.0f + exp((actualVm + 5.0f) / 5.0f));
        real gamma_d = 1.0f / (1.0f + exp((50.0f - actualVm) / 20.0f));

        real f_inf = 1.0f / (1.0f + exp((actualVm + 20.0f) / 7.0f));
        real alpha_f = 1102.5f * exp(-(actualVm + 27.0f) * (actualVm + 27.0f) / 225.0f);
        real beta_f = 200.0f / (1.0f + exp((13.0f - actualVm) / 10.0f));
        real gamma_f = 180.0f / (1.0f + exp((actualVm + 30.0f) / 10.0f)) + 20.0f;

        real f2_inf = 0.67f / (1.0f + exp((actualVm + 35.0f) / 7.0f)) + 0.33f;
        real alpha_f2; // !!!
        alpha_f2 = 562.0f * exp(-(actualVm + 27.0f) * (actualVm + 27.0f) / 240.0f);
        real beta_f2 = 31.0f / (1.0f + exp((25.0f - actualVm) / 10.0f));
        real gamma_f2; // !!!
        gamma_f2 = 80.0f / (1.0f + exp((30.0f + actualVm) / 10.0f));

        real fCaSS_inf = 0.6f / (1.0f + (actualCa_SS * actualCa_SS * 400.0f)) + 0.4f;
        real tau_fCaSS = 80.0f / (1.0f + (actualCa_SS * actualCa_SS * 400.0f)) + 2.0f;

#if defined(EPI) || defined(MCELL)

        real s_inf = 1.0f / (1.0f + exp((actualVm + 20.0f) / 5.0f));
        real tau_s = 85.0f * exp(-(actualVm + 45.0f) * (actualVm + 45.0f) / 320.0f) + 5.0f / (1.0f + exp((actualVm - 20.0f) / 5.0f)) + 3.0f;

#endif // EPI || MCELL
#ifdef ENDO

        real s_inf = 1.0f / (1.0f + exp((actualVm + 28.0f) / 5.0f));
        real tau_s = 1000.0f * exp(-(actualVm + 67.0f) * (actualVm + 67.0f) / 1000.0f) + 8.0f;

#endif // ENDO

        real r_inf = 1.0f / (1.0f + exp((20.0f - actualVm) / 6.0f));
        real tau_r = 9.5f * exp(-(actualVm + 40.0f) * (actualVm + 40.0f) / 1800.0f) + 0.8f;

        // Rush-Larsen method - update approximations
        d_X_r1[index] = X_r1_inf - (X_r1_inf - actualX_r1) * exp(-delta_t / (alpha_X_r1 * beta_X_r1));
        d_X_r2[index] = X_r2_inf - (X_r2_inf - actualX_r2) * exp(-delta_t / (alpha_X_r2 * beta_X_r2));
        d_X_s[index] = X_s_inf - (X_s_inf - actualX_s) * exp(-delta_t / (alpha_X_s * beta_X_s + 80.0f));
        d_m[index] = m_inf - (m_inf - actualm) * exp(-delta_t / (alpha_m * beta_m));
        d_h[index] = h_inf - (h_inf - actualh) * exp(-delta_t * (alpha_h + beta_h));
        d_j[index] = j_inf - (j_inf - actualj) * exp(-delta_t * (alpha_j + beta_j));
        d_d[index] = inf - (inf - actuald) * exp(-delta_t / (alpha_d * beta_d + gamma_d));
        d_f[index] = f_inf - (f_inf - actualf) * exp(-delta_t / (alpha_f + beta_f + gamma_f));
        d_f2[index] = f2_inf - (f2_inf - actualf2) * exp(-delta_t / (alpha_f2 + beta_f2 + gamma_f2));
        d_fCaSS[index] = fCaSS_inf - (fCaSS_inf - actualfCaSS) * exp(-delta_t / tau_fCaSS);
        d_s[index] = s_inf - (s_inf - actuals) * exp(-delta_t / tau_s);
        d_r[index] = r_inf - (r_inf - actualr) * exp(-delta_t / tau_r);

        // Explicit method - auxiliary variables
        real Ileak = V_leak * (actualCa_SR - actualCa_i);
        real Iup = V_maxup / (1.0f + (K_up * K_up) / (actualCa_i * actualCa_i));
        real k_CaSR = max_SR - (max_SR - min_SR) / (1.0f + (EC / actualCa_SR) * (EC / actualCa_SR));
        real k1 = k1_prime / k_CaSR;
        real k2 = k2_prime * k_CaSR;
        real O = k1 * actualCa_SS * actualCa_SS * actualR_prime / (k3 + k1 * actualCa_SS * actualCa_SS);
        real Irel = V_rel * O * (actualCa_SR - actualCa_SS);
        real Ixfer = V_xfer * (actualCa_SS - actualCa_i);
        real Ca_i_bufC; // !!!
        Ca_i_bufC = 1.0f / (1.0f + bufC * K_bufC / ((K_bufC + actualCa_i) * (K_bufC + actualCa_i)));
        real Ca_SR_bufSR; // !!!
        Ca_SR_bufSR = 1.0f / (1.0f + bufSR * K_bufSR / ((K_bufSR + actualCa_SR) * (K_bufSR + actualCa_SR)));
        real Ca_SS_bufSS; // !!!
        Ca_SS_bufSS = 1.0f / (1.0f + bufSS * K_bufSS / ((K_bufSS + actualCa_SS) * (K_bufSS + actualCa_SS)));

        // Explicit method - RHS of the state variables
        real RHS_R_prime_term = (k4 * (1.0f - actualR_prime)) - k2 * actualCa_SS * actualR_prime;
        real RHS_Ca_i_term = Ca_i_bufC * ((((Ileak - Iup) * V_SR / V_C) + Ixfer) - (((IbCa + IpCa) - 2.0f * INaCa) * Cm / (2.0f * V_C * F)));
        real RHS_Ca_SR_term = Ca_SR_bufSR * (Iup - Ileak - Irel);
        real RHS_Ca_SS_term = Ca_SS_bufSS * (((-ICaL * Cm / (2.0 * V_SS * F)) + (Irel * V_SR / V_SS)) - (Ixfer * V_C / V_SS));
        real RHS_Na_i_term = -((INa + IbNa + 3.0f * INaK + 3.0f * INaCa) * Cm / (V_C * F));
        real RHS_K_i_term = -((IK1 + Ito + IKr + IKs - 2.0f * INaK + IpK + stim) * Cm / (V_C * F));

        // Explicit method - update approximations
        d_R_prime[index] = actualR_prime + (delta_t * RHS_R_prime_term);
        d_Ca_i[index] = actualCa_i + (delta_t * RHS_Ca_i_term);
        d_Ca_SR[index] = actualCa_SR + (delta_t * RHS_Ca_SR_term);
        d_Ca_SS[index] = actualCa_SS + (delta_t * RHS_Ca_SS_term);
        d_Na_i[index] = actualNa_i + (delta_t * RHS_Na_i_term);
        d_K_i[index] = actualK_i + (delta_t * RHS_K_i_term);

#endif // FE
    }
}

#endif // TT2

#ifdef MV

__global__ void computeApprox(int Nx, int Ny, real delta_t, real phi_x, real phi_y, real diff_coeff, real actualTime, Stimulus *d_stimuli, real *d_Vm, real *d_partRHS, real *d_v, real *d_w, real *d_s)
{
    // Obtain the index of the thread
    int i = blockIdx.y * blockDim.y + threadIdx.y; // Coordinate y
    int j = blockIdx.x * blockDim.x + threadIdx.x; // Coordinate x

    if (i < Ny && j < Nx)
    {
        // Calculate the index in the 1D array
        int index = i * Nx + j;

        // Stimulation
        real stim = 0.0f;

#pragma unroll
        for (int si = 0; si < numberOfStimuli; si++)
        {
            const Stimulus &stimulus = d_stimuli[si];
            if (actualTime >= stimulus.begin && actualTime <= stimulus.begin + stimulus.duration &&
                j >= stimulus.xMinDisc && j <= stimulus.xMaxDisc &&
                i >= stimulus.yMinDisc && i <= stimulus.yMaxDisc)
            {
                stim = stimulus.amplitude;
                break;
            }
        }

        // Get the actual central value
        real actualVm = d_Vm[index];

        // State variable
        real actualv = d_v[index];
        real actualw = d_w[index];
        real actuals = d_s[index];

#if defined(SSIADI) || defined(THETASSIADI) || defined(FE)

        // Boundary conditions
        real im1Vm = (i - 1 == -1) ? d_Vm[index + Nx] : d_Vm[index - Nx];
        real ip1Vm = (i + 1 == Ny) ? d_Vm[index - Nx] : d_Vm[index + Nx];
        real jm1Vm = (j - 1 == -1) ? d_Vm[index + 1] : d_Vm[index - 1];
        real jp1Vm = (j + 1 == Nx) ? d_Vm[index - 1] : d_Vm[index + 1];

        real diff_term = diff_coeff * (phi_x * (jm1Vm - 2.0f * actualVm + jp1Vm) + phi_y * (im1Vm - 2.0f * actualVm + ip1Vm));

#endif // SSIADI || THETASSIADI || FE

        // Calculate RHS of the equations
        // Auxiliary variables
        real Htheta_w = (actualVm - theta_w > 0.0f) ? 1.0f : 0.0f;
        real Htheta_o = (actualVm - theta_o > 0.0f) ? 1.0f : 0.0f;
        real Htheta_v = (actualVm - theta_v > 0.0f) ? 1.0f : 0.0f;
        real Htheta_vminus = (actualVm - theta_vminus > 0.0f) ? 1.0f : 0.0f;

        real tau_o = (1.0f - Htheta_o) * tau_o1 + Htheta_o * tau_o2;
        real tau_so = tau_so1 + (tau_so2 - tau_so1) * (1.0f + tanh(k_so * (actualVm - u_so))) * 0.5f;
        real tau_vminus = (1.0f - Htheta_vminus) * tau_v1minus + Htheta_vminus * tau_v2minus;

        // Currents
        real J_fi = -actualv * Htheta_v * (actualVm - theta_v) * (u_u - actualVm) / tau_fi;
        real J_so = ((actualVm - u_o) * (1.0f - Htheta_w) / tau_o) + (Htheta_w / tau_so);
        real J_si = -Htheta_w * actualw * actuals / tau_si;

        // RHS of the state variables
        real RHS_Vm_term = J_fi + J_so + J_si;
        real RHS_v_term = (1.0f - Htheta_v) * (((actualVm < theta_vminus) ? 1.0f : 0.0f) - actualv) / tau_vminus - (Htheta_v * actualv / tau_vplus);
        real RHS_w_term = (1.0f - Htheta_w) * (((1.0f - Htheta_o) * (1.0f - (actualVm / tau_winf)) + Htheta_o * w_infstar) - actualw) / (tau_w1minus + (((tau_w2minus - tau_w1minus) * (1.0f + tanh(k_wminus * (actualVm - u_wminus)))) * 0.5f)) - (Htheta_w * actualw / tau_wplus);
        real RHS_s_term = (((1.0f + tanh(k_s * (actualVm - u_s))) * 0.5f) - actuals) / ((1.0f - Htheta_w) * tau_s1 + Htheta_w * tau_s2);

#if defined(SSIADI) || defined(THETASSIADI)

        // Calculate Vmtilde -> utilde = u^n + 0.5 * dt * (A*u^n + R(u^n))
        real Vmtilde = actualVm + 0.5f * diff_term + (0.5f * delta_t * (stim - RHS_Vm_term));

        // Auxiliary variables for Rush-Larsen or Forward Euler
        real v_inf = ((actualVm < theta_vminus) ? 1.0f : 0.0f);
        real tau_v_RL = (tau_vplus * tau_vminus) / (tau_vplus - tau_vplus * Htheta_v + tau_vminus * Htheta_v);
        real v_inf_RL = (tau_vplus * v_inf * (1.0f - Htheta_v)) / (tau_vplus - tau_vplus * Htheta_v + tau_vminus * Htheta_v);

        real w_inf = ((1.0f - Htheta_o) * (1.0f - (actualVm / tau_winf)) + Htheta_o * w_infstar);
        real tau_wminus = tau_w1minus + (tau_w2minus - tau_w1minus) * (1.0f + tanh(k_wminus * (actualVm - u_wminus))) * 0.5f;
        real tau_w_RL = (tau_wplus * tau_wminus) / (tau_wplus - tau_wplus * Htheta_w + tau_wminus * Htheta_w);
        real w_inf_RL = (tau_wplus * w_inf * (1.0f - Htheta_w)) / (tau_wplus - tau_wplus * Htheta_w + tau_wminus * Htheta_w);

        real tau_s = (1.0f - Htheta_w) * tau_s1 + Htheta_w * tau_s2;
        real s_inf_RL = (1.0f + tanh(k_s * (actualVm - u_s))) * 0.5f;

        // Calculate approximations with Rush-Larsen or Forward Euler using half time step
        real vtilde, wtilde, stilde;
        (tau_v_RL > 1.0e-10)
            ? (vtilde = v_inf_RL - (v_inf_RL - actualv) * exp(-0.5f * delta_t / tau_v_RL))
            : (vtilde = actualv + 0.5f * delta_t * (1.0f - Htheta_v) * (v_inf - actualv) / tau_vminus - Htheta_v * actualv / tau_vplus);

        (tau_w_RL > 1.0e-10)
            ? (wtilde = w_inf_RL - (w_inf_RL - actualw) * exp(-0.5f * delta_t / tau_w_RL))
            : (wtilde = actualw + 0.5f * delta_t * (1.0f - Htheta_w) * (w_inf - actualw) / tau_wminus - Htheta_w * actualw / tau_wplus);

        (tau_s > 1.0e-10)
            ? (stilde = s_inf_RL - (s_inf_RL - actuals) * exp(-0.5f * delta_t / tau_s))
            : (stilde = actuals + 0.5f * delta_t * (s_inf_RL - actuals) / tau_s);

        // Calculate RHS of the equations with approximations
        // Auxiliary variables
        Htheta_w = (Vmtilde - theta_w > 0.0f) ? 1.0f : 0.0f;
        Htheta_o = (Vmtilde - theta_o > 0.0f) ? 1.0f : 0.0f;
        Htheta_v = (Vmtilde - theta_v > 0.0f) ? 1.0f : 0.0f;
        Htheta_vminus = (Vmtilde - theta_vminus > 0.0f) ? 1.0f : 0.0f;

        tau_o = (1.0f - Htheta_o) * tau_o1 + Htheta_o * tau_o2;
        tau_so = tau_so1 + (tau_so2 - tau_so1) * (1.0f + tanh(k_so * (Vmtilde - u_so))) * 0.5f;
        tau_vminus = (1.0f - Htheta_vminus) * tau_v1minus + Htheta_vminus * tau_v2minus;

        // Currents
        real J_fi_tilde = -actualv * Htheta_v * (Vmtilde - theta_v) * (u_u - Vmtilde) / tau_fi;
        real J_so_tilde = ((Vmtilde - u_o) * (1.0f - Htheta_w) / ((1.0f - Htheta_o) * tau_o1 + Htheta_o * tau_o2)) + (Htheta_w / (tau_so1 + (((tau_so2 - tau_so1) * (1.0f + tanh(k_so * (Vmtilde - u_so)))) * 0.5f)));
        real J_si_tilde = -Htheta_w * wtilde * stilde / tau_si;

        // Update d_partRHS
        real RHS_Vmtilde_term = J_fi_tilde + J_so_tilde + J_si_tilde;
        d_partRHS[index] = delta_t * (stim - RHS_Vmtilde_term);

        // Update auxiliary variables for Rush-Larsen or Forward Euler with approximations
        v_inf = ((Vmtilde < theta_vminus) ? 1.0f : 0.0f);
        tau_v_RL = (tau_vplus * tau_vminus) / (tau_vplus - tau_vplus * Htheta_v + tau_vminus * Htheta_v);
        v_inf_RL = (tau_vplus * v_inf * (1.0f - Htheta_v)) / (tau_vplus - tau_vplus * Htheta_v + tau_vminus * Htheta_v);

        w_inf = ((1.0f - Htheta_o) * (1.0f - (Vmtilde / tau_winf)) + Htheta_o * w_infstar);
        tau_wminus = tau_w1minus + (tau_w2minus - tau_w1minus) * (1.0f + tanh(k_wminus * (Vmtilde - u_wminus))) * 0.5f;
        tau_w_RL = (tau_wplus * tau_wminus) / (tau_wplus - tau_wplus * Htheta_w + tau_wminus * Htheta_w);
        w_inf_RL = (tau_wplus * w_inf * (1.0f - Htheta_w)) / (tau_wplus - tau_wplus * Htheta_w + tau_wminus * Htheta_w);

        tau_s = (1.0f - Htheta_w) * tau_s1 + Htheta_w * tau_s2;
        s_inf_RL = (1.0f + tanh(k_s * (Vmtilde - u_s))) * 0.5f;

        // Update state variables with Rush-Larsen or RK2 (Heun's method) using approximations
        (tau_v_RL > 1.0e-10)
            ? (d_v[index] = v_inf_RL - (v_inf_RL - actualv) * exp(-delta_t / tau_v_RL))
            : (d_v[index] = actualv + delta_t * (1.0f - Htheta_v) * (v_inf - vtilde) / tau_vminus - Htheta_v * vtilde / tau_vplus);

        (tau_w_RL > 1.0e-10)
            ? (d_w[index] = w_inf_RL - (w_inf_RL - actualw) * exp(-delta_t / tau_w_RL))
            : (d_w[index] = actualw + delta_t * (1.0f - Htheta_w) * (w_inf - wtilde) / tau_wminus - Htheta_w * wtilde / tau_wplus);

        (tau_s > 1.0e-10)
            ? (d_s[index] = s_inf_RL - (s_inf_RL - actuals) * exp(-delta_t / tau_s))
            : (d_s[index] = actuals + delta_t * (s_inf_RL - stilde) / tau_s);

#endif // SSIADI || THETASSIADI

#ifdef OSADI

        // Update d_partRHS
        d_partRHS[index] = delta_t * (stim - RHS_Vm_term);

        // Auxiliary variables for Rush-Larsen or Forward Euler
        real v_inf = ((actualVm < theta_vminus) ? 1.0f : 0.0f);
        real tau_v_RL = (tau_vplus * tau_vminus) / (tau_vplus - tau_vplus * Htheta_v + tau_vminus * Htheta_v);
        real v_inf_RL = (tau_vplus * v_inf * (1.0f - Htheta_v)) / (tau_vplus - tau_vplus * Htheta_v + tau_vminus * Htheta_v);

        real w_inf = ((1.0f - Htheta_o) * (1.0f - (actualVm / tau_winf)) + Htheta_o * w_infstar);
        real tau_wminus = tau_w1minus + (tau_w2minus - tau_w1minus) * (1.0f + tanh(k_wminus * (actualVm - u_wminus))) * 0.5f;
        real tau_w_RL = (tau_wplus * tau_wminus) / (tau_wplus - tau_wplus * Htheta_w + tau_wminus * Htheta_w);
        real w_inf_RL = (tau_wplus * w_inf * (1.0f - Htheta_w)) / (tau_wplus - tau_wplus * Htheta_w + tau_wminus * Htheta_w);

        real tau_s = (1.0f - Htheta_w) * tau_s1 + Htheta_w * tau_s2;
        real s_inf_RL = (1.0f + tanh(k_s * (actualVm - u_s))) * 0.5f;

        // Update state variables with Rush-Larsen or Forward Euler
        (tau_v_RL > 1.0e-10)
            ? (d_v[index] = v_inf_RL - (v_inf_RL - actualv) * exp(-delta_t / tau_v_RL))
            : (d_v[index] = actualv + delta_t * (1.0f - Htheta_v) * (v_inf - actualv) / tau_vminus - Htheta_v * actualv / tau_vplus);

        (tau_w_RL > 1.0e-10)
            ? (d_w[index] = w_inf_RL - (w_inf_RL - actualw) * exp(-delta_t / tau_w_RL))
            : (d_w[index] = actualw + delta_t * (1.0f - Htheta_w) * (w_inf - actualw) / tau_wminus - Htheta_w * actualw / tau_wplus);

        (tau_s > 1.0e-10)
            ? (d_s[index] = s_inf_RL - (s_inf_RL - actuals) * exp(-delta_t / tau_s))
            : (d_s[index] = actuals + delta_t * (s_inf_RL - actuals) / tau_s);

#endif // OSADI

#ifdef FE

        // Update d_partRHS (auxiliary variable) with Forward Euler
        d_partRHS[index] = actualVm + diff_term + delta_t * (stim - RHS_Vm_term);

        // Auxiliary variables for Rush-Larsen or Forward Euler
        real v_inf = ((actualVm < theta_vminus) ? 1.0f : 0.0f);
        real tau_v_RL = (tau_vplus * tau_vminus) / (tau_vplus - tau_vplus * Htheta_v + tau_vminus * Htheta_v);
        real v_inf_RL = (tau_vplus * v_inf * (1.0f - Htheta_v)) / (tau_vplus - tau_vplus * Htheta_v + tau_vminus * Htheta_v);

        real w_inf = ((1.0f - Htheta_o) * (1.0f - (actualVm / tau_winf)) + Htheta_o * w_infstar);
        real tau_wminus = tau_w1minus + (tau_w2minus - tau_w1minus) * (1.0f + tanh(k_wminus * (actualVm - u_wminus))) * 0.5f;
        real tau_w_RL = (tau_wplus * tau_wminus) / (tau_wplus - tau_wplus * Htheta_w + tau_wminus * Htheta_w);
        real w_inf_RL = (tau_wplus * w_inf * (1.0f - Htheta_w)) / (tau_wplus - tau_wplus * Htheta_w + tau_wminus * Htheta_w);

        real tau_s = (1.0f - Htheta_w) * tau_s1 + Htheta_w * tau_s2;
        real s_inf_RL = (1.0f + tanh(k_s * (actualVm - u_s))) * 0.5f;

        // Update state variables with Rush-Larsen or Forward Euler
        (tau_v_RL > 1.0e-10)
            ? (d_v[index] = v_inf_RL + (actualv - v_inf_RL) * exp(-delta_t / tau_v_RL))
            : (d_v[index] = actualv + delta_t * ((1.0f - Htheta_v) * (v_inf - actualv) / tau_vminus - Htheta_v * actualv / tau_vplus));

        (tau_w_RL > 1.0e-10)
            ? (d_w[index] = w_inf_RL + (actualw - w_inf_RL) * exp(-delta_t / tau_w_RL))
            : (d_w[index] = actualw + delta_t * ((1.0f - Htheta_w) * (w_inf - actualw) / tau_wminus - Htheta_w * actualw / tau_wplus));

        (tau_s > 1.0e-10)
            ? (d_s[index] = s_inf_RL + (actuals - s_inf_RL) * exp(-delta_t / tau_s))
            : (d_s[index] = actuals + delta_t * (s_inf_RL - actuals) / tau_s);

#endif // FE

#ifdef HV // Hundsdorfer-Verwer scheme (for anisotropy)

        // Boundary conditions
        real im1Vm = (i - 1 == -1) ? d_Vm[index + Nx] : d_Vm[index - Nx];
        real ip1Vm = (i + 1 == Ny) ? d_Vm[index - Nx] : d_Vm[index + Nx];
        real jm1Vm = (j - 1 == -1) ? d_Vm[index + 1] : d_Vm[index - 1];
        real jp1Vm = (j + 1 == Nx) ? d_Vm[index - 1] : d_Vm[index + 1];
        // TODO: solve cross derivatives in boundaries

        real diff_term = diff_coeff * (phi_x * (jm1Vm - 2.0f * actualVm + jp1Vm) + phi_y * (im1Vm - 2.0f * actualVm + ip1Vm));

#endif // HV
    }
}

#endif // MV

#if defined(SSIADI) || defined(THETASSIADI) || defined(OSADI)

__global__ void prepareRHSiDiff(int Nx, int Ny, real phi_y, real diff_coeff, real tau, real *d_Vm, real *d_RHS, real *d_partRHS)
{
    // Obtain the index of the thread
    int i = blockIdx.y * blockDim.y + threadIdx.y; // Coordinate y
    int j = blockIdx.x * blockDim.x + threadIdx.x; // Coordinate x

    if (i < Ny && j < Nx)
    {
        // Calculate the index in the 1D array
        int index = i * Nx + j;

        // Load central
        real actualVm = d_Vm[index];

        // Boundary conditions
        real im1Vm = (i - 1 == -1) ? d_Vm[index + Nx] : d_Vm[index - Nx];
        real ip1Vm = (i + 1 == Ny) ? d_Vm[index - Nx] : d_Vm[index + Nx];

        // Calculate the diffusion term
        real diff_term = diff_coeff * tau * phi_y * (im1Vm - 2.0f * actualVm + ip1Vm);

        // Update d_RHS
        d_RHS[index] = actualVm + diff_term + 0.5f * d_partRHS[index]; // this 0.5f is associated to a two dimension case of ADI
    }
}

__global__ void prepareRHSjDiff(int Nx, int Ny, real phi_x, real diff_coeff, real tau, real *d_Vm, real *d_RHS, real *d_partRHS)
{
    // Obtain the index of the thread
    int i = blockIdx.y * blockDim.y + threadIdx.y; // Coordinate y
    int j = blockIdx.x * blockDim.x + threadIdx.x; // Coordinate x

    if (i < Ny && j < Nx)
    {
        // Calculate the index in the 1D array
        int index = i * Nx + j;

        // Load central
        real actualVm = d_Vm[index];

        // Boundary conditions
        real jm1Vm = (j - 1 == -1) ? d_Vm[index + 1] : d_Vm[index - 1];
        real jp1Vm = (j + 1 == Nx) ? d_Vm[index - 1] : d_Vm[index + 1];

        // Calculate the diffusion term
        real diff_term = diff_coeff * tau * phi_x * (jm1Vm - 2.0f * actualVm + jp1Vm);

        // Update d_RHS
        d_RHS[index] = actualVm + diff_term + 0.5f * d_partRHS[index]; // this 0.5f is associated to a two dimension case of ADI
    }
}

#endif // SSIADI || THETASSIADI

#ifdef OSADI

__global__ void prepareRHS(int Nx, int Ny, real *d_Vm, real *d_partRHS)
{
    // Obtain the index of the thread
    int i = blockIdx.y * blockDim.y + threadIdx.y; // Coordinate y
    int j = blockIdx.x * blockDim.x + threadIdx.x; // Coordinate x

    if (i < Ny && j < Nx)
    {
        int index = i * Nx + j;
        d_Vm[index] = d_Vm[index] + 0.5f * d_partRHS[index]; // this 0.5f is associated to a two dimension case of ADI
    }
}

#endif // OSADI

__global__ void parallelThomasVertical(int numSys, int sysSize, real *d, real *la, real *lb, real *lc)
{
    // Each thread will solve a system
    int systemIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (systemIdx < numSys)
    {
        // Local variables
        real c_prime[MAX_SYS_SIZE];
        real d_prime[MAX_SYS_SIZE];

        // Declare and load shared memory
        // __shared__ real shared_la[MAX_SYS_SIZE];
        // __shared__ real shared_lb[MAX_SYS_SIZE];
        // __shared__ real shared_lc[MAX_SYS_SIZE];

        // for (int i = 0; i < sysSize && threadIdx.x == 0; i++)
        // {
        //     shared_la[i] = la[i];
        //     shared_lb[i] = lb[i];
        //     shared_lc[i] = lc[i];
        // }
        // __syncthreads();

        c_prime[0] = lc[0] / lb[0];
        d_prime[0] = d[systemIdx] / lb[0];

        for (int i = 1; i < sysSize; i++)
        {
            real denom = 1.0f / (lb[i] - c_prime[i - 1] * la[i]);
            if (i < sysSize - 1)
            {
                c_prime[i] = lc[i] * denom;
            }
            d_prime[i] = (d[systemIdx + i * numSys] - d_prime[i - 1] * la[i]) * denom;
        }

        d[systemIdx + (sysSize - 1) * numSys] = d_prime[sysSize - 1];
        for (int i = sysSize - 2; i >= 0; i--)
        {
            d[systemIdx + i * numSys] = d_prime[i] - c_prime[i] * d[systemIdx + (i + 1) * numSys];
        }
    }
}

__global__ void parallelThomasHorizontal(int numSys, int sysSize, real *d, real *la, real *lb, real *lc)
{
    // Each thread will solve a system
    int systemIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (systemIdx < numSys)
    {
        // Calculate the offset of the system
        int offset = systemIdx * sysSize;

        // Local variables
        real c_prime[MAX_SYS_SIZE];
        real d_prime[MAX_SYS_SIZE];

        // Declare and load shared memory
        // __shared__ real shared_la[MAX_SYS_SIZE];
        // __shared__ real shared_lb[MAX_SYS_SIZE];
        // __shared__ real shared_lc[MAX_SYS_SIZE];

        // for (int i = 0; i < sysSize && threadIdx.x == 0; i++)
        // {
        //     shared_la[i] = la[i];
        //     shared_lb[i] = lb[i];
        //     shared_lc[i] = lc[i];
        // }
        // __syncthreads();

        c_prime[0] = lc[0] / lb[0];
        d_prime[0] = d[offset] / lb[0];

        for (int i = 1; i < sysSize; i++)
        {
            real denom = 1.0f / (lb[i] - c_prime[i - 1] * la[i]);
            if (i < sysSize - 1)
            {
                c_prime[i] = lc[i] * denom;
            }
            d_prime[i] = (d[offset + i] - d_prime[i - 1] * la[i]) * denom;
        }

        d[offset + sysSize - 1] = d_prime[sysSize - 1];
        for (int i = sysSize - 2; i >= 0; i--)
        {
            d[offset + i] = d_prime[i] - c_prime[i] * d[offset + i + 1];
        }
    }
}

// TODO: Implement the following functions that use prefactorization
// __global__ void parallelThomasVertical(int numSys, int sysSize, real *d, real *la, real *c_prime, real *denominator)
// {
//     // Each thread will solve a system
//     int systemIdx = blockIdx.x * blockDim.x + threadIdx.x;

//     if (systemIdx < numSys)
//     {
//         // Local variables
//         real d_prime[MAX_SYS_SIZE];

//         // Declare and load shared memory
//         __shared__ real shared_la[MAX_SYS_SIZE];
//         __shared__ real shared_c_prime[MAX_SYS_SIZE];
//         __shared__ real shared_denominator[MAX_SYS_SIZE];

//         for (int i = 0; i < sysSize && threadIdx.x == 0; i++)
//         {
//             shared_la[i] = la[i];
//             shared_c_prime[i] = c_prime[i];
//             shared_denominator[i] = denominator[i];
//         }
//         __syncthreads();

//         d_prime[0] = d[systemIdx] * shared_denominator[0];

//         for (int i = 1; i < sysSize; i++)
//         {
//             d_prime[i] = (d[systemIdx + i * numSys] - d_prime[i - 1] * shared_la[i]) * shared_denominator[i];
//         }

//         d[systemIdx + (sysSize - 1) * numSys] = d_prime[sysSize - 1];
//         for (int i = sysSize - 2; i >= 0; i--)
//         {
//             d[systemIdx + i * numSys] = d_prime[i] - shared_c_prime[i] * d[systemIdx + (i + 1) * numSys];
//         }
//     }
// }

// __global__ void parallelThomasHorizontal(int numSys, int sysSize, real *d, real *la, real *c_prime, real *denominator)
// {
//     // Each thread will solve a system
//     int systemIdx = blockIdx.x * blockDim.x + threadIdx.x;

//     if (systemIdx < numSys)
//     {
//         // Calculate the offset of the system
//         int offset = systemIdx * sysSize;

//         // Local variables
//         real d_prime[MAX_SYS_SIZE];

//         // Declare and load shared memory
//         __shared__ real shared_la[MAX_SYS_SIZE];
//         __shared__ real shared_c_prime[MAX_SYS_SIZE];
//         __shared__ real shared_denominator[MAX_SYS_SIZE];

//         for (int i = 0; i < sysSize && threadIdx.x == 0; i++)
//         {
//             shared_la[i] = la[i];
//             shared_c_prime[i] = c_prime[i];
//             shared_denominator[i] = denominator[i];
//         }
//         __syncthreads();

//         d_prime[0] = d[offset] * shared_denominator[0];

//         for (int i = 1; i < sysSize; i++)
//         {
//             d_prime[i] = (d[offset + i] - d_prime[i - 1] * shared_la[i]) * shared_denominator[i];
//         }

//         d[offset + sysSize - 1] = d_prime[sysSize - 1];
//         for (int i = sysSize - 2; i >= 0; i--)
//         {
//             d[offset + i] = d_prime[i] - shared_c_prime[i] * d[offset + i + 1];
//         }
//     }
// }

#endif // USE_CUDA

#endif // GPU_FUNCTIONS_H