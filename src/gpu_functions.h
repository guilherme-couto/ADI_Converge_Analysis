#ifndef GPU_FUNCTIONS_H
#define GPU_FUNCTIONS_H

#include "simulation_config.h"

#ifdef MONODOMAIN

#ifdef AFHN

// Kernel to compute the approximate solution of the reaction-diffusion system and update the state variables
__global__ void computeApprox(int Nx, int Ny, real delta_t, real phi_x, real phi_y, real diff_coeff, real actualTime, real *d_V, real *d_partRHS, real *d_W, Stimulus *d_stimuli)
{
    // Obtain the index of the thread
    int i = blockIdx.y * blockDim.y + threadIdx.y; // Coordinate y
    int j = blockIdx.x * blockDim.x + threadIdx.x; // Coordinate x

    if (i < Ny && j < Nx)
    {
        // Calculate the index in the 1D array
        int index = i * Nx + j;
        real actualV = d_V[index];

        // State variable
        real actualW = d_W[index];

        // Calculate RHS_V at actual time
        real RHS_V_term = ((G * actualV * (1.0f - (actualV / vth)) * (1.0f - (actualV / vp))) + (eta1 * actualV * actualW)) / (Cm * chi);

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
                stim = stimulus.strength;
                break;
            }
        }

#if defined(SSIADI) || defined(THETASSIADI)

        int index_im1 = (i - 1) * Nx + j;
        int index_ip1 = (i + 1) * Nx + j;
        int index_jm1 = i * Nx + (j - 1);
        int index_jp1 = i * Nx + (j + 1);

        // Boundary conditions
        (i - 1 == -1) ? (index_im1 = index_ip1) : index_im1;
        (i + 1 == Ny) ? (index_ip1 = index_im1) : index_ip1;
        (j - 1 == -1) ? (index_jm1 = index_jp1) : index_jm1;
        (j + 1 == Nx) ? (index_jp1 = index_jm1) : index_jp1;

        real im1V = d_V[index_im1];
        real ip1V = d_V[index_ip1];
        real jm1V = d_V[index_jm1];
        real jp1V = d_V[index_jp1];

        real diff_term = diff_coeff * (phi_x * (jm1V - 2.0f * actualV + jp1V) + phi_y * (im1V - 2.0f * actualV + ip1V));

        // Calculate Vtilde -> utilde = u^n + 0.5 * dt * (A*u^n + R(u^n))
        real Vtilde = actualV + 0.5f * diff_term + (0.5f * delta_t * (stim - RHS_V_term));

        // Calculate approximation for state variables
        real RHS_W_term = eta2 * ((actualV / vp) - (eta3 * actualW));
        real Wtilde = actualW + (0.5f * delta_t * RHS_W_term);

        // Preparing part of the RHS of the following linear systems
        real RHS_Vtilde_term = ((G * Vtilde * (1.0f - (Vtilde / vth)) * (1.0f - (Vtilde / vp))) + (eta1 * Vtilde * Wtilde)) / (Cm * chi); // RHS with V* and W*
        d_partRHS[index] = delta_t * (stim - RHS_Vtilde_term);

        // Update state variables with RK2 -> Wn+1 = Wn + dt*R(V*, W*)
        real RHS_Wtilde_term = eta2 * ((Vtilde / vp) - (eta3 * Wtilde));
        d_W[index] = actualW + delta_t * RHS_Wtilde_term;

#endif // SSIADI || THETASSIADI

#ifdef OSADI

        // Calculate part of the RHS of the following linear systems with Forward Euler
        d_partRHS[index] = delta_t * (stim - RHS_V_term);

        // Update state variables
        d_W[index] = actualW + delta_t * eta2 * ((actualV / vp) - (eta3 * actualW)); // with Forward Euler -> Wn+1 = Wn + dt*R(Vn, Wn)

#endif // OSADI

#ifdef FE

        // Update variables explicitly
        real diff_term = diff_coeff * (phi_x * (d_V[lim(i, j - 1, Nx)] - 2.0f * actualV + d_V[lim(i, j + 1, Nx)]) + phi_y * (d_V[lim(i - 1, j, Ny)] - 2.0f * actualV + d_V[lim(i + 1, j, Ny)]));
        d_partRHS[index] = actualV + diff_term + delta_t * (stim - RHS_V_term);
        d_W[index] = actualW + delta_t * eta2 * ((actualV / vp) - (eta3 * actualW));

#endif // FE
        
    }
}

#endif // AFHN

#ifdef TT2

// Kernel to compute the approximate solution of the reaction-diffusion system and update the state variables
__global__ void computeApprox(int Nx, int Ny, real delta_t, real phi_x, real phi_y, real diff_coeff, real actualTime, real *d_V, real *d_partRHS, real *d_X_r1, real *d_X_r2, real *d_X_s, real *d_m, real *d_h, real *d_j, real *d_d, real *d_f, real *d_f2, real *d_fCaSS, real *d_s, real *d_r, real *d_Ca_i, real *d_Ca_SR, real *d_Ca_SS, real *d_R_prime, real *d_Na_i, real *d_K_i, Stimulus *d_stimuli)
{
    // Obtain the index of the thread
    int i = blockIdx.y * blockDim.y + threadIdx.y; // Coordinate y
    int j = blockIdx.x * blockDim.x + threadIdx.x; // Coordinate x

    if (i < Ny && j < Nx)
    {
        // Calculate the index in the 1D array
        int index = i * Nx + j;

        // Boundary conditions
        int index_im1 = (i - 1) * Nx + j;
        int index_ip1 = (i + 1) * Nx + j;
        int index_jm1 = i * Nx + (j - 1);
        int index_jp1 = i * Nx + (j + 1);

        (i - 1 == -1) ? (index_im1 = index_ip1) : index_im1;
        (i + 1 == Ny) ? (index_ip1 = index_im1) : index_ip1;
        (j - 1 == -1) ? (index_jm1 = index_jp1) : index_jm1;
        (j + 1 == Nx) ? (index_jp1 = index_jm1) : index_jp1;

        real actualV = d_V[index];
        real im1V = d_V[index_im1];
        real ip1V = d_V[index_ip1];
        real jm1V = d_V[index_jm1];
        real jp1V = d_V[index_jp1];

        real diff_term = diff_coeff * (phi_x * (jm1V - 2.0f * actualV + jp1V) + phi_y * (im1V - 2.0f * actualV + ip1V));

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
        real VmENa = actualV - (RTONF * log(Na_o / actualNa_i));
        real E_K = (RTONF * log(K_o / actualK_i));
        real VmEK = actualV - E_K;
        real alpha_K1 = 0.1f / (1.0f + exp(0.06f * (actualV - E_K - 200.0f)));
        real beta_K1 = (3.0f * exp(0.0002f * (actualV - E_K + 100.0f)) + exp(0.1f * (actualV - E_K - 10.0f))) / (1.0f + exp(-0.5f * (actualV - E_K)));
        real E_Ks = (RTONF * log((K_o + p_KNa * Na_o) / (actualK_i + p_KNa * actualNa_i)));
        real E_Ca = 0.5f * RTONF * log(Ca_o / actualCa_i);

        // Currents
        real INa = G_Na * (actualm * actualm * actualm) * actualh * actualj * VmENa;
        real IbNa = G_bNa * VmENa;
        real IK1 = G_K1 * (alpha_K1 / (alpha_K1 + beta_K1)) * VmEK;
        real Ito = G_to * actualr * actuals * VmEK;
        real IKr = G_Kr * sqrt(K_o / 5.4f) * actualX_r1 * actualX_r2 * VmEK;
        real IKs = G_Ks * actualX_s * actualX_s * (actualV - E_Ks);
        real ICaL; // !!!
        (actualV < 15.0f - 1.0e-5f)
            ? (ICaL = G_CaL * actuald * actualf * actualf2 * actualfCaSS * 4.0f * (actualV - 15.0f) * (F * F) * (0.25f * actualCa_SS * exp(2.0f * (actualV - 15.0f) * FONRT) - Ca_o) / (R * T * (exp(2.0f * (actualV - 15.0f) * FONRT) - 1.0f)))
            : (ICaL = G_CaL * actuald * actualf * actualf2 * actualfCaSS * 2.0f * F * (0.25f * actualCa_SS - Ca_o));
        real INaK = ((((p_KNa * K_o) / (K_o + K_mK)) * actualNa_i) / (actualNa_i + K_mNa)) / (1.0f + (0.1245f * exp(((-0.1f) * actualV * FONRT))) + (0.0353f * exp(((-actualV) * FONRT))));
        real INaCa; // !!!
        INaCa = (k_NaCa * ((exp((gamma_I_NaCa * actualV * FONRT)) * (actualNa_i * actualNa_i * actualNa_i) * Ca_o) - (exp(((gamma_I_NaCa - 1.0f) * actualV * FONRT)) * (Na_o * Na_o * Na_o) * actualCa_i * alpha))) / (((K_mNa_i * K_mNa_i * K_mNa_i) + (Na_o * Na_o * Na_o)) * (K_mCa + Ca_o) * (1.0f + (k_sat * exp(((gamma_I_NaCa)*actualV * FONRT)))));
        real IpCa = (G_pCa * actualCa_i) / (K_pCa + actualCa_i);
        real IpK = (G_pK * VmEK) / (1.0f + exp((25.0f - actualV) / 5.98f));
        real IbCa = G_bCa * (actualV - E_Ca);

        // RHS_V at actual time
        real RHS_V_term = INa + IbNa + IK1 + Ito + IKr + IKs + ICaL + INaK + INaCa + IpCa + IpK + IbCa;

        real stim = 0.0f;
        for (int si = 0; si < numberOfStimuli; si++)
        {
            if (actualTime >= d_stimuli[si].begin && actualTime <= d_stimuli[si].begin + d_stimuli[si].duration && j >= d_stimuli[si].xMinDisc && j <= d_stimuli[si].xMaxDisc && i >= d_stimuli[si].yMinDisc && i <= d_stimuli[si].yMaxDisc)
            {
                stim = d_stimuli[si].strength;
                break;
            }
        }

        // Calculate Vtilde -> utilde = u^n + 0.5 * dt * (A*u^n + R(u^n))
        real Vtilde = actualV + 0.5f * diff_term + (0.5f * delta_t * (stim - RHS_V_term));

        // Preparing part of the RHS of the following linear systems
        // Calculate approximation for state variables
        // Rush-Larsen method - auxiliary variables
        real X_r1_inf = 1.0f / (1.0f + exp((-26.0f - actualV) / 7.0f));
        real alpha_X_r1 = 450.0f / (1.0f + exp((-45.0f - actualV) / 10.0f));
        real beta_X_r1 = 6.0f / (1.0f + exp((30.0f + actualV) / 11.5f));

        real X_r2_inf = 1.0f / (1.0f + exp((actualV + 88.0f) / 24.0f));
        real alpha_X_r2 = 3.0f / (1.0f + exp((-60.0f - actualV) / 20.0f));
        real beta_X_r2 = 1.12f / (1.0f + exp((actualV - 60.0f) / 20.0f));

        real X_s_inf = 1.0f / (1.0f + exp((-5.0f - actualV) / 14.0f));
        real alpha_X_s = 1400.0f / sqrt(1.0f + exp((5.0f - actualV) / 6.0f));
        real beta_X_s = 1.0f / (1.0f + exp((-35.0f + actualV) / 15.0f));

        real m_inf = 1.0f / ((1.0f + exp((-56.86 - actualV) / 9.03f)) * (1.0f + exp((-56.86f - actualV) / 9.03f)));
        real alpha_m = 1.0f / (1.0f + exp((-60.0f - actualV) / 5.0f));
        real beta_m = 0.1f / (1.0f + exp((actualV + 35.0f) / 5.0f)) + (0.1f / (1.0f + exp((actualV - 50.0f) / 200.0f)));

        real h_inf = 1.0f / ((1.0f + exp((actualV + 71.55f) / 7.43f)) * (1.0f + exp((actualV + 71.55f) / 7.43f)));
        real alpha_h;
        (actualV < -40.0f)
            ? (alpha_h = 0.057f * exp(-(actualV + 80.0f) / 6.8f))
            : (alpha_h = 0.0f);
        real beta_h;
        (actualV < -40.0f)
            ? (beta_h = 2.7f * exp(0.079f * actualV) + 3.1f * 1.0e5f * exp(0.3485f * actualV))
            : (beta_h = 0.77f / (0.13f * (1.0f + exp((actualV + 10.66f) / -11.1f))));

        real j_inf = 1.0f / ((1.0f + exp((actualV + 71.55f) / 7.43f)) * (1.0f + exp((actualV + 71.55f) / 7.43f)));
        real alpha_j;
        (actualV < -40.0f)
            ? (alpha_j = ((-25428.0f * exp(0.2444f * actualV) - (6.948e-6f * exp((-0.04391f) * actualV))) * (actualV + 37.78f)) / (1.0f + exp(0.311f * (actualV + 79.23f))))
            : (alpha_j = 0.0f);
        real beta_j;
        (actualV < -40.0f)
            ? (beta_j = (0.02424f * exp(-0.01052f * actualV)) / (1.0f + exp(-0.1378f * (actualV + 40.14f))))
            : (beta_j = (0.6f * exp(0.057f * actualV)) / (1.0f + exp(-0.1f * (actualV + 32.0f))));

        real inf = 1.0f / (1.0f + exp((-8.0f - actualV) / 7.5f));
        real alpha_d = 1.4f / (1.0f + exp((-35.0f - actualV) / 13.0f)) + 0.25f;
        real beta_d = 1.4f / (1.0f + exp((actualV + 5.0f) / 5.0f));
        real gamma_d = 1.0f / (1.0f + exp((50.0f - actualV) / 20.0f));

        real f_inf = 1.0f / (1.0f + exp((actualV + 20.0f) / 7.0f));
        real alpha_f = 1102.5f * exp(-(actualV + 27.0f) * (actualV + 27.0f) / 225.0f);
        real beta_f = 200.0f / (1.0f + exp((13.0f - actualV) / 10.0f));
        real gamma_f = 180.0f / (1.0f + exp((actualV + 30.0f) / 10.0f)) + 20.0f;

        real f2_inf = 0.67f / (1.0f + exp((actualV + 35.0f) / 7.0f)) + 0.33f;
        real alpha_f2; // !!!
        alpha_f2 = 562.0f * exp(-(actualV + 27.0f) * (actualV + 27.0f) / 240.0f);
        real beta_f2 = 31.0f / (1.0f + exp((25.0f - actualV) / 10.0f));
        real gamma_f2; // !!!
        gamma_f2 = 80.0f / (1.0f + exp((30.0f + actualV) / 10.0f));

        real fCaSS_inf = 0.6f / (1.0f + (actualCa_SS * actualCa_SS * 400.0f)) + 0.4f;
        real tau_fCaSS = 80.0f / (1.0f + (actualCa_SS * actualCa_SS * 400.0f)) + 2.0f;

#if defined(EPI) || defined(MCELL)

        real s_inf = 1.0f / (1.0f + exp((actualV + 20.0f) / 5.0f));
        real tau_s = 85.0f * exp(-(actualV + 45.0f) * (actualV + 45.0f) / 320.0f) + 5.0f / (1.0f + exp((actualV - 20.0f) / 5.0f)) + 3.0f;

#endif // EPI || MCELL
#ifdef ENDO

        real s_inf = 1.0f / (1.0f + exp((actualV + 28.0f) / 5.0f));
        real tau_s = 1000.0f * exp(-(actualV + 67.0f) * (actualV + 67.0f) / 1000.0f) + 8.0f;

#endif // ENDO

        real r_inf = 1.0f / (1.0f + exp((20.0f - actualV) / 6.0f));
        real tau_r = 9.5f * exp(-(actualV + 40.0f) * (actualV + 40.0f) / 1800.0f) + 0.8f;

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

        // Explicit method - update approximations
        real R_primetilde = actualR_prime + (0.5f * delta_t * RHS_R_prime_term);
        real Ca_itilde = actualCa_i + (0.5f * delta_t * RHS_Ca_i_term);
        real Ca_SRtilde = actualCa_SR + (0.5f * delta_t * RHS_Ca_SR_term);
        real Ca_SStilde = actualCa_SS + (0.5f * delta_t * RHS_Ca_SS_term);
        real Na_itilde = actualNa_i + (0.5f * delta_t * RHS_Na_i_term);
        real K_itilde = actualK_i + (0.5f * delta_t * RHS_K_i_term);

        // Auxiliary variables with Vtilde
        real VmENatilde = Vtilde - (RTONF * log(Na_o / Na_itilde));
        real E_Ktilde = (RTONF * log(K_o / K_itilde));
        real VmEKtilde = Vtilde - E_Ktilde;
        real alpha_K1tilde = 0.1f / (1.0f + exp(0.06f * (Vtilde - E_Ktilde - 200.0f)));
        real beta_K1tilde = (3.0f * exp(0.0002f * (Vtilde - E_Ktilde + 100.0f)) + exp(0.1f * (Vtilde - E_Ktilde - 10.0f))) / (1.0f + exp(-0.5f * (Vtilde - E_Ktilde)));
        real E_Kstilde = (RTONF * log((K_o + p_KNa * Na_o) / (K_itilde + p_KNa * Na_itilde)));
        real E_Catilde = 0.5f * RTONF * log(Ca_o / Ca_itilde);

        // Currents with Vtilde
        real INatilde = G_Na * (mtilde * mtilde * mtilde) * htilde * jtilde * VmENatilde;
        real IbNatilde = G_bNa * VmENatilde;
        real IK1tilde = G_K1 * (alpha_K1tilde / (alpha_K1tilde + beta_K1tilde)) * VmEKtilde;
        real Itotilde = G_to * rtilde * stilde * VmEKtilde;
        real IKrtilde = G_Kr * sqrt(K_o / 5.4f) * X_r1tilde * X_r2tilde * VmEKtilde;
        real IKstilde = G_Ks * X_stilde * X_stilde * (Vtilde - E_Kstilde);
        real ICaLtilde; // !!!
        (Vtilde < 15.0f - 1.0e-5f)
            ? (ICaLtilde = G_CaL * dtilde * ftilde * f2tilde * fCaSStilde * 4.0f * (Vtilde - 15.0f) * (F * F) * (0.25f * Ca_SStilde * exp(2.0f * (Vtilde - 15.0f) * FONRT) - Ca_o) / (R * T * (exp(2.0f * (Vtilde - 15.0f) * FONRT) - 1.0f)))
            : (ICaLtilde = G_CaL * dtilde * ftilde * f2tilde * fCaSStilde * 2.0f * F * (0.25f * Ca_SStilde - Ca_o));
        real INaKtilde = ((((p_KNa * K_o) / (K_o + K_mK)) * Na_itilde) / (Na_itilde + K_mNa)) / (1.0f + (0.1245f * exp(((-0.1f) * Vtilde * FONRT))) + (0.0353f * exp(((-Vtilde) * FONRT))));
        real INaCatilde; // !!!
        INaCatilde = (k_NaCa * ((exp((gamma_I_NaCa * Vtilde * FONRT)) * (Na_itilde * Na_itilde * Na_itilde) * Ca_o) - (exp(((gamma_I_NaCa - 1.0f) * Vtilde * FONRT)) * (Na_o * Na_o * Na_o) * Ca_itilde * alpha))) / (((K_mNa_i * K_mNa_i * K_mNa_i) + (Na_o * Na_o * Na_o)) * (K_mCa + Ca_o) * (1.0f + (k_sat * exp(((gamma_I_NaCa)*Vtilde * FONRT)))));
        real IpCatilde = (G_pCa * Ca_itilde) / (K_pCa + Ca_itilde);
        real IpKtilde = (G_pK * VmEKtilde) / (1.0f + exp((25.0f - Vtilde) / 5.98f));
        real IbCatilde = G_bCa * (Vtilde - E_Catilde);

        // RHS of the main equation with Vtilde
        real RHS_Vtilde_term = INatilde + IbNatilde + IK1tilde + Itotilde + IKrtilde + IKstilde + ICaLtilde + INaKtilde + INaCatilde + IpCatilde + IpKtilde + IbCatilde;
        d_partRHS[index] = delta_t * (stim - RHS_Vtilde_term);

        // Update state variables
        // RHS of the state variables with tilde approximations
        // Rush-Larsen method - auxiliary variables
        real X_r1_inftilde = 1.0f / (1.0f + exp((-26.0f - Vtilde) / 7.0f));
        real alpha_X_r1tilde = 450.0f / (1.0f + exp((-45.0f - Vtilde) / 10.0f));
        real beta_X_r1tilde = 6.0f / (1.0f + exp((30.0f + Vtilde) / 11.5f));

        real X_r2_inftilde = 1.0f / (1.0f + exp((Vtilde + 88.0f) / 24.0f));
        real alpha_X_r2tilde = 3.0f / (1.0f + exp((-60.0f - Vtilde) / 20.0f));
        real beta_X_r2tilde = 1.12f / (1.0f + exp((Vtilde - 60.0f) / 20.0f));

        real X_s_inftilde = 1.0f / (1.0f + exp((-5.0f - Vtilde) / 14.0f));
        real alpha_X_stilde = 1400.0f / sqrt(1.0f + exp((5.0f - Vtilde) / 6.0f));
        real beta_X_stilde = 1.0f / (1.0f + exp((-35.0f + Vtilde) / 15.0f));

        real m_inftilde = 1.0f / ((1.0f + exp((-56.86 - Vtilde) / 9.03f)) * (1.0f + exp((-56.86f - Vtilde) / 9.03f)));
        real alpha_mtilde = 1.0f / (1.0f + exp((-60.0f - Vtilde) / 5.0f));
        real beta_mtilde = 0.1f / (1.0f + exp((Vtilde + 35.0f) / 5.0f)) + (0.1f / (1.0f + exp((Vtilde - 50.0f) / 200.0f)));

        real h_inftilde = 1.0f / ((1.0f + exp((Vtilde + 71.55f) / 7.43f)) * (1.0f + exp((Vtilde + 71.55f) / 7.43f)));
        real alpha_htilde;
        (Vtilde < -40.0f)
            ? (alpha_htilde = 0.057f * exp(-(Vtilde + 80.0f) / 6.8f))
            : (alpha_htilde = 0.0f);
        real beta_htilde;
        (Vtilde < -40.0f)
            ? (beta_htilde = 2.7f * exp(0.079f * Vtilde) + 3.1f * 1.0e5f * exp(0.3485f * Vtilde))
            : (beta_htilde = 0.77f / (0.13f * (1.0f + exp((Vtilde + 10.66f) / -11.1f))));

        real j_inftilde = 1.0f / ((1.0f + exp((Vtilde + 71.55f) / 7.43f)) * (1.0f + exp((Vtilde + 71.55f) / 7.43f)));
        real alpha_jtilde;
        (Vtilde < -40.0f)
            ? (alpha_jtilde = ((-25428.0f * exp(0.2444f * Vtilde) - (6.948e-6f * exp((-0.04391f) * Vtilde))) * (Vtilde + 37.78f)) / (1.0f + exp(0.311f * (Vtilde + 79.23f))))
            : (alpha_jtilde = 0.0f);
        real beta_jtilde;
        (Vtilde < -40.0f)
            ? (beta_jtilde = (0.02424f * exp(-0.01052f * Vtilde)) / (1.0f + exp(-0.1378f * (Vtilde + 40.14f))))
            : (beta_jtilde = (0.6f * exp(0.057f * Vtilde)) / (1.0f + exp(-0.1f * (Vtilde + 32.0f))));

        real d_inftilde = 1.0f / (1.0f + exp((-8.0f - Vtilde) / 7.5f));
        real alpha_dtilde = 1.4f / (1.0f + exp((-35.0f - Vtilde) / 13.0f)) + 0.25f;
        real beta_dtilde = 1.4f / (1.0f + exp((Vtilde + 5.0f) / 5.0f));
        real gamma_dtilde = 1.0f / (1.0f + exp((50.0f - Vtilde) / 20.0f));

        real f_inftilde = 1.0f / (1.0f + exp((Vtilde + 20.0f) / 7.0f));
        real alpha_ftilde = 1102.5f * exp(-(Vtilde + 27.0f) * (Vtilde + 27.0f) / 225.0f);
        real beta_ftilde = 200.0f / (1.0f + exp((13.0f - Vtilde) / 10.0f));
        real gamma_ftilde = 180.0f / (1.0f + exp((Vtilde + 30.0f) / 10.0f)) + 20.0f;

        real f2_inftilde = 0.67f / (1.0f + exp((Vtilde + 35.0f) / 7.0f)) + 0.33f;
        real alpha_f2tilde; // !!!
        alpha_f2tilde = 562.0f * exp(-(Vtilde + 27.0f) * (Vtilde + 27.0f) / 240.0f);
        real beta_f2tilde = 31.0f / (1.0f + exp((25.0f - Vtilde) / 10.0f));
        real gamma_f2tilde; // !!!
        gamma_f2tilde = 80.0f / (1.0f + exp((30.0f + Vtilde) / 10.0f));

        real fCaSS_inftilde = 0.6f / (1.0f + (Ca_SStilde * Ca_SStilde * 400.0f)) + 0.4f;
        real tau_fCaSStilde = 80.0f / (1.0f + (Ca_SStilde * Ca_SStilde * 400.0f)) + 2.0f;

#if defined(EPI) || defined(MCELL)

        real s_inftilde = 1.0f / (1.0f + exp((Vtilde + 20.0f) / 5.0f));
        real tau_stilde = 85.0f * exp(-(Vtilde + 45.0f) * (Vtilde + 45.0f) / 320.0f) + 5.0f / (1.0f + exp((Vtilde - 20.0f) / 5.0f)) + 3.0f;

#endif // EPI || MCELL
#ifdef ENDO

        real s_inftilde = 1.0f / (1.0f + exp((Vtilde + 28.0f) / 5.0f));
        real tau_stilde = 1000.0f * exp(-(Vtilde + 67.0f) * (Vtilde + 67.0f) / 1000.0f) + 8.0f;

#endif // ENDO

        real r_inftilde = 1.0f / (1.0f + exp((20.0f - Vtilde) / 6.0f));
        real tau_rtilde = 9.5f * exp(-(Vtilde + 40.0f) * (Vtilde + 40.0f) / 1800.0f) + 0.8f;

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

        // Explicit method - update approximations
        d_R_prime[index] = actualR_prime + (delta_t * RHS_R_primetilde_term);
        d_Ca_i[index] = actualCa_i + (delta_t * RHS_Ca_itilde_term);
        d_Ca_SR[index] = actualCa_SR + (delta_t * RHS_Ca_SRtilde_term);
        d_Ca_SS[index] = actualCa_SS + (delta_t * RHS_Ca_SStilde_term);
        d_Na_i[index] = actualNa_i + (delta_t * RHS_Na_itilde_term);
        d_K_i[index] = actualK_i + (delta_t * RHS_K_itilde_term);
    }
}

#endif // TT2

#if defined(SSIADI) || defined(THETASSIADI)

__global__ void prepareRHSiDiff(int Nx, int Ny, real phi_y, real diff_coeff, real tau, real *d_V, real *d_RHS, real *d_partRHS)
{
    // Obtain the index of the thread
    int i = blockIdx.y * blockDim.y + threadIdx.y; // Coordinate y
    int j = blockIdx.x * blockDim.x + threadIdx.x; // Coordinate x

    if (i < Ny && j < Nx)
    {
        // Calculate the index in the 1D array
        int index = i * Nx + j;

        int index_im1 = (i - 1) * Nx + j;
        int index_ip1 = (i + 1) * Nx + j;

        (i - 1 == -1) ? (index_im1 = index_ip1) : index_im1;
        (i + 1 == Ny) ? (index_ip1 = index_im1) : index_ip1;

        real actualV = d_V[index];

        real diff_term = diff_coeff * tau * phi_y * (d_V[index_im1] - 2.0f * actualV + d_V[index_ip1]);
        d_RHS[index] = actualV + diff_term + 0.5f * d_partRHS[index]; // this 0.5f is associated to a two dimension case of ADI
    }
}

__global__ void prepareRHSjDiff(int Nx, int Ny, real phi_x, real diff_coeff, real tau, real *d_V, real *d_RHS, real *d_partRHS)
{
    // Obtain the index of the thread
    int i = blockIdx.y * blockDim.y + threadIdx.y; // Coordinate y
    int j = blockIdx.x * blockDim.x + threadIdx.x; // Coordinate x

    if (i < Ny && j < Nx)
    {
        // Calculate the index in the 1D array
        int index = i * Nx + j;

        int index_jm1 = i * Nx + (j - 1);
        int index_jp1 = i * Nx + (j + 1);

        (j - 1 == -1) ? (index_jm1 = index_jp1) : index_jm1;
        (j + 1 == Nx) ? (index_jp1 = index_jm1) : index_jp1;

        real actualV = d_V[index];

        real diff_term = diff_coeff * tau * phi_x * (d_V[index_jm1] - 2.0f * actualV + d_V[index_jp1]);
        d_RHS[index] = actualV + diff_term + 0.5f * d_partRHS[index]; // this 0.5f is associated to a two dimension case of ADI
    }
}

#endif // SSIADI || THETASSIADI

#ifdef OSADI

__global__ void prepareRHS(int Nx, int Ny, real *d_V, real *d_partRHS)
{
    // Obtain the index of the thread
    int i = blockIdx.y * blockDim.y + threadIdx.y; // Coordinate y
    int j = blockIdx.x * blockDim.x + threadIdx.x; // Coordinate x

    if (i < Ny && j < Nx)
    {
        int index = i * Nx + j;
        d_V[index] = d_V[index] + 0.5f * d_partRHS[index]; // this 0.5f is associated to a two dimension case of ADI
    }
}

#endif // OSADI

// From GLOSTER, Andrew et al. Efficient Interleaved Batch Matrix Solvers for CUDA. arXiv preprint arXiv:1909.04539, 2019.
// Kernel to solve the tridiagonal system using the Thomas algorithm
// numSys -> Number of systems to solve
// sysSize -> Size of each system
// d -> Result vector
// la -> Lower diagonal
// lb -> Main diagonal
// lc -> Upper diagonal
__global__ void parallelThomasVertical(int numSys, int sysSize, real *d, real *la, real *lb, real *lc)
{
    int previousRow, nextRow;
    int currentRow = blockIdx.x * blockDim.x + threadIdx.x;
    int i = 0;

    if (currentRow < numSys)
    {
        // 1st: update auxiliary arrays
        d[currentRow] = d[currentRow] / lb[i];

#pragma unroll
        for (i = 1; i < sysSize; i++)
        {
            previousRow = currentRow;
            currentRow += numSys;

            d[currentRow] = (d[currentRow] - la[i] * d[previousRow]) / (lb[i]);
        }

        // 2nd: update solution
        d[currentRow] = d[currentRow];

#pragma unroll
        for (i = sysSize - 2; i >= 0; i--)
        {
            nextRow = currentRow;
            currentRow -= numSys;

            d[currentRow] = d[currentRow] - lc[i] * d[nextRow];
        }
    }
}

__global__ void parallelThomasHorizontal(int numSys, int sysSize, real *d, real *la, real *lb, real *lc)
{
    int previousColumn, nextColumn;
    int currentColumn = blockIdx.x * blockDim.x + threadIdx.x;
    currentColumn *= numSys;
    int i = 0;

    if (currentColumn < numSys)
    {
        // 1st: update auxiliary arrays
        d[currentColumn] = d[currentColumn] / lb[i];

        for (i = 1; i < sysSize; i++)
        {
            previousColumn = currentColumn;
            currentColumn += 1;

            d[currentColumn] = (d[currentColumn] - la[i] * d[previousColumn]) / (lb[i]);
        }

        // 2nd: update solution
        d[currentColumn] = d[currentColumn];

        for (i = sysSize - 2; i >= 0; i--)
        {
            nextColumn = currentColumn;
            currentColumn -= 1;

            d[currentColumn] = d[currentColumn] - lc[i] * d[nextColumn];
        }
    }
}

#endif // MONODOMAIN

#endif // GPU_FUNCTIONS_H