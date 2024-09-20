#ifndef CPU_METHODS_H
#define CPU_METHODS_H

#include "auxfuncs.h"

void runSimulation(char *method, real delta_t, real delta_x, real theta)
{
    // Number of steps
    int N = round(L / delta_x) + 1;     // Spatial steps (square tissue)
    int M = round(totalTime / delta_t); // Number of time steps

    // Allocate and populate time array
    real *time = (real *)malloc(M * sizeof(real));
    initializeTimeArray(time, M, delta_t);
#ifndef CABLEEQ
    // Allocate 2D arrays for variables
    real **V, **Vtilde, **RHS, **partRHS, **exact;
    V = (real **)malloc(N * sizeof(real *));
    Vtilde = (real **)malloc(N * sizeof(real *));
    RHS = (real **)malloc(N * sizeof(real *));
    partRHS = (real **)malloc(N * sizeof(real *));
    exact = (real **)malloc(N * sizeof(real *));
#else
    // Allocate 1D arrays for variables
    real *V, *Vtilde, *RHS, *partRHS, *exact;
    V = (real *)malloc(N * sizeof(real));
    Vtilde = (real *)malloc(N * sizeof(real));
    RHS = (real *)malloc(N * sizeof(real));
    partRHS = (real *)malloc(N * sizeof(real));
    exact = (real *)malloc(N * sizeof(real));
#endif // not CABLEEQ
    real *c_prime = (real *)malloc(N * sizeof(real)); // aux for Thomas
    real *d_prime = (real *)malloc(N * sizeof(real)); // aux for Thomas
    real *d = (real *)malloc(N * sizeof(real));
    real *result = (real *)malloc(N * sizeof(real));
#ifdef MONODOMAIN
#ifdef AFHN
    real **W = (real **)malloc(N * sizeof(real *));
#endif // AFHN
#endif // MONODOMAIN
#ifdef CABLEEQ
#ifdef AFHN
    real *W = (real *)malloc(N * sizeof(real));
#endif // AFHN
#endif // CABLEEQ
#ifndef CABLEEQ
    for (int i = 0; i < N; i++)
    {
        V[i] = (real *)malloc(N * sizeof(real));
        Vtilde[i] = (real *)malloc(N * sizeof(real));
        RHS[i] = (real *)malloc(N * sizeof(real));
        partRHS[i] = (real *)malloc(N * sizeof(real));
        exact[i] = (real *)malloc(N * sizeof(real));
#ifdef MONODOMAIN
#ifdef AFHN
        W[i] = (real *)malloc(N * sizeof(real));
#endif // AFHN
#endif // MONODOMAIN
    }
#endif // not CABLEEQ

#ifdef CONVERGENCE_ANALYSIS
    initialize2DVariableWithExactSolution(V, N, delta_x);
#else
#ifndef CABLEEQ
    initialize2DVariableWithValue(V, N, V_init);
#else
    initialize1DVariableWithValue(V, N, V_init);
#endif // not CABLEEQ
#endif // CONVERGENCE_ANALYSIS

#ifdef MONODOMAIN
#ifdef AFHN
    initialize2DVariableWithValue(W, N, 0.0f);
#endif // AFHN
#endif // MONODOMAIN
#ifdef CABLEEQ
#ifdef AFHN
    initialize1DVariableWithValue(W, N, 0.0f);
#endif // AFHN
#endif // CABLEEQ

#ifdef INIT_WITH_SPIRAL
    char *pathToSpiralFiles = (char *)malloc(MAX_STRING_SIZE * sizeof(char));
    snprintf(pathToSpiralFiles, MAX_STRING_SIZE * sizeof(char), "./spiral_files/%s/%s/%s/lastV_0.0005_0.0005.txt", REAL_TYPE, PROBLEM, CELL_MODEL);
    initialize2DVariableFromFile(V, N, pathToSpiralFiles, delta_x, "V");
    snprintf(pathToSpiralFiles, MAX_STRING_SIZE * sizeof(char), "./spiral_files/%s/%s/%s/lastW_0.0005_0.0005.txt", REAL_TYPE, PROBLEM, CELL_MODEL);
    initialize2DVariableFromFile(W, N, pathToSpiralFiles, delta_x, "W");
    free(pathToSpiralFiles);
#endif // INIT_WITH_SPIRAL

    // Auxiliary arrays for Thomas algorithm
    real *la = (real *)malloc(N * sizeof(real)); // subdiagonal
    real *lb = (real *)malloc(N * sizeof(real)); // diagonal
    real *lc = (real *)malloc(N * sizeof(real)); // superdiagonal

    // Populate auxiliary arrays for Thomas algorithm
    real phi = (delta_t / (delta_x * delta_x));
    if (strcmp(method, "ADI") == 0 || strcmp(method, "SSI-ADI") == 0)
    {
        populateDiagonalThomasAlgorithm(la, lb, lc, N, 0.5f * phi);
    }
    else if (strcmp(method, "theta-ADI") == 0)
    {
        populateDiagonalThomasAlgorithm(la, lb, lc, N, theta * phi);
    }

    // Create directories
    char *pathToSaveData = (char *)malloc(MAX_STRING_SIZE * sizeof(char));
    createDirectories(method, theta, pathToSaveData);

    int timeStepCounter = 0;
    real actualTime = 0.0f;

    // Measure total execution time
    real startTime = omp_get_wtime();
#ifndef CABLEEQ
    if (strcmp(method, "ADI") == 0)
    {
        while (timeStepCounter < M)
        {
            // Get time step
            actualTime = time[timeStepCounter];

            // ================================================!
            //  Calcula V em n + 1/2 -> Resultado vai para RHS !
            // ================================================!
            real x, y;
            for (int j = 0; j < N; j++)
            {
                for (int i = 0; i < N; i++)
                {
                    x = j * delta_x;
                    y = i * delta_x;
#if defined(LINMONO) || defined(DIFF)
                    d[i] = 0.5f * phi * V[i][lim(j - 1, N)] + (1.0f - 2.0f * 0.5f * phi) * V[i][j] + 0.5f * phi * V[i][lim(j + 1, N)] + 0.5f * delta_t * forcingTerm(x, y, actualTime);
#endif // LINMONO || DIFF
                }

                tridiag(la, lb, lc, c_prime, d_prime, N, d, result);
                for (int i = 0; i < N; i++)
                {
                    RHS[i][j] = result[i];
                }
            }

            // ================================================!
            //  Calcula V em n + 1 -> Resultado vai para V     !
            // ================================================!
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    x = j * delta_x;
                    y = i * delta_x;
#if defined(LINMONO) || defined(DIFF)
                    d[j] = 0.5f * phi * RHS[lim(i - 1, N)][j] + (1.0f - 2.0f * 0.5f * phi) * RHS[i][j] + 0.5f * phi * RHS[lim(i + 1, N)][j] + 0.5f * delta_t * forcingTerm(x, y, actualTime + delta_t);
#endif // LINMONO || DIFF
                }

                tridiag(la, lb, lc, c_prime, d_prime, N, d, result);
                for (int j = 0; j < N; j++)
                {
                    V[i][j] = result[j];
                }
            }

            // Update time step counter
            timeStepCounter++;
        }
    }

#if defined(LINMONO) || defined(MONODOMAIN)
    else if (strcmp(method, "SSI-ADI") == 0)
    {
        while (timeStepCounter < M)
        {
            // Get time step
            actualTime = time[timeStepCounter];

            // ================================================!
            //  Calcula Approx.                                !
            // ================================================!
            real x, y;
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    x = j * delta_x;
                    y = i * delta_x;

#ifdef LINMONO
                    real diff_term = (sigma / (chi * Cm)) * 0.5f * phi * (V[i][lim(j - 1, N)] + V[lim(i - 1, N)][j] - 4.0f * V[i][j] + V[i][lim(j + 1, N)] + V[lim(i + 1, N)][j]);
                    real for_term = forcingTerm(x, y, actualTime + (0.5f * delta_t)) / (chi * Cm);
                    real reac_term = G * V[i][j] / Cm;
                    Vtilde[i][j] = V[i][j] + diff_term + (0.5f * delta_t * (for_term - reac_term));

                    // Preparing part of the RHS of the following linear systems
                    real reac_tilde_term = G * Vtilde[i][j] / Cm;
                    partRHS[i][j] = delta_t * (for_term - reac_tilde_term);
#endif // LINMONO

#ifdef MONODOMAIN
                    real diff_term = (sigma / (chi * Cm)) * 0.5f * phi * (V[i][lim(j - 1, N)] + V[lim(i - 1, N)][j] - 4.0f * V[i][j] + V[i][lim(j + 1, N)] + V[lim(i + 1, N)][j]);

#ifdef AFHN
                    real for_term = forcingTerm(x, y, actualTime + (0.5f * delta_t), W[i][j]) / (chi * Cm);
                    real RHS_V_term = RHS_V(V[i][j], W[i][j]);
#endif // AFHN

                    Vtilde[i][j] = V[i][j] + diff_term + (0.5f * delta_t * (for_term - RHS_V_term));

// Preparing part of the RHS of the following linear systems
#ifdef AFHN
                    // Calculate W approximation
                    // real RHS_W_term = RHS_W(V[i][j], W[i][j]);
                    // real Wtilde = W[i][j] + (0.5f * delta_t * RHS_W_term);
                    real RHS_Vtilde_term = RHS_V(Vtilde[i][j], W[i][j]);
#endif // AFHN

                    partRHS[i][j] = delta_t * (for_term - RHS_Vtilde_term);

#ifdef AFHN
                    // Update Wn+1
                    real RHS_W_term = RHS_W(V[i][j], W[i][j]);
                    real Wtilde = W[i][j] + (0.5f * delta_t * RHS_W_term);
                    W[i][j] = W[i][j] + delta_t * RHS_W(Vtilde[i][j], Wtilde);
                    // RHS_Vtilde_term = RHS_V(Vtilde[i][j], Wtilde);
#endif // AFHN
                    // partRHS[i][j] = delta_t * (for_term - RHS_Vtilde_term);
#endif // MONODOMAIN
                }
            }

            // ================================================!
            //  Calcula V em n + 1/2 -> Resultado vai para RHS !
            // ================================================!
            for (int j = 0; j < N; j++)
            {
                for (int i = 0; i < N; i++)
                {
                    real diff_term = (sigma / (chi * Cm)) * 0.5f * phi * (V[i][lim(j - 1, N)] - 2.0f * V[i][j] + V[i][lim(j + 1, N)]);
                    d[i] = V[i][j] + diff_term + 0.5f * partRHS[i][j];
                }

                tridiag(la, lb, lc, c_prime, d_prime, N, d, result);
                for (int i = 0; i < N; i++)
                {
                    RHS[i][j] = result[i];
                }
            }

            // ================================================!
            //  Calcula V em n + 1 -> Resultado vai para V     !
            // ================================================!
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    real diff_term = (sigma / (chi * Cm)) * 0.5f * phi * (RHS[lim(i - 1, N)][j] - 2.0f * RHS[i][j] + RHS[lim(i + 1, N)][j]);
                    d[j] = RHS[i][j] + diff_term + 0.5f * partRHS[i][j];
                }

                tridiag(la, lb, lc, c_prime, d_prime, N, d, result);
                for (int j = 0; j < N; j++)
                {
                    V[i][j] = result[j];
                }
            }

            // Update time step counter
            timeStepCounter++;
        }
    }

    else if (strcmp(method, "theta-ADI") == 0)
    {
        while (timeStepCounter < M)
        {
            // Get time step
            actualTime = time[timeStepCounter];

            // ================================================!
            //  Calcula Approx.                                !
            // ================================================!
            real x, y;
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    x = j * delta_x;
                    y = i * delta_x;

#ifdef LINMONO
                    real diff_term = (sigma / (chi * Cm)) * phi * (V[i][lim(j - 1, N)] + V[lim(i - 1, N)][j] - 4.0f * V[i][j] + V[i][lim(j + 1, N)] + V[lim(i + 1, N)][j]);
                    real for_term = forcingTerm(x, y, actualTime + (0.5f * delta_t)) / (chi * Cm);
                    real reac_term = G * V[i][j] / Cm;
                    Vtilde[i][j] = V[i][j] + diff_term + (delta_t * (for_term - reac_term));

                    // Preparing part of the RHS of the following linear systems
                    real reac_tilde_term = G * Vtilde[i][j] / Cm;
                    partRHS[i][j] = delta_t * (for_term - ((1.0f - theta) * reac_term) - (theta * reac_tilde_term));
#endif // LINMONO

#ifdef MONODOMAIN
                    real diff_term = (sigma / (chi * Cm)) * phi * (V[i][lim(j - 1, N)] + V[lim(i - 1, N)][j] - 4.0f * V[i][j] + V[i][lim(j + 1, N)] + V[lim(i + 1, N)][j]);

#ifdef AFHN
                    real for_term = forcingTerm(x, y, actualTime + (0.5f * delta_t), W[i][j]) / (chi * Cm);
                    real RHS_V_term = RHS_V(V[i][j], W[i][j]);
#endif // AFHN

                    Vtilde[i][j] = V[i][j] + diff_term + (delta_t * (for_term - RHS_V_term));

// Preparing part of the RHS of the following linear systems
#ifdef AFHN
                    real RHS_Vtilde_term = RHS_V(Vtilde[i][j], W[i][j]);
#endif // AFHN

                    partRHS[i][j] = delta_t * (for_term - ((1.0f - theta) * RHS_V_term) - (theta * RHS_Vtilde_term));

#ifdef AFHN
                    // Update Wn+1
                    real RHS_W_term = RHS_W(V[i][j], W[i][j]);
                    real Wtilde = W[i][j] + (delta_t * RHS_W_term);
                    W[i][j] = W[i][j] + delta_t * RHS_W(Vtilde[i][j], Wtilde);
                    // RHS_Vtilde_term = RHS_V(Vtilde[i][j], Wtilde);
#endif // AFHN
                    // partRHS[i][j] = delta_t * (for_term - ((1.0f - theta) * RHS_V_term) - (theta * RHS_Vtilde_term));
#endif // MONODOMAIN
                }
            }

            // ================================================!
            //  Calcula V em n + 1/2 -> Resultado vai para RHS !
            // ================================================!
            for (int j = 0; j < N; j++)
            {
                for (int i = 0; i < N; i++)
                {
                    real diff_term = (sigma / (chi * Cm)) * (1.0f - theta) * phi * (V[i][lim(j - 1, N)] - 2.0f * V[i][j] + V[i][lim(j + 1, N)]);
                    d[i] = V[i][j] + diff_term + 0.5f * partRHS[i][j];
                }

                tridiag(la, lb, lc, c_prime, d_prime, N, d, result);
                for (int i = 0; i < N; i++)
                {
                    RHS[i][j] = result[i];
                }
            }

            // ================================================!
            //  Calcula V em n + 1 -> Resultado vai para V     !
            // ================================================!
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    real diff_term = (sigma / (chi * Cm)) * (1.0f - theta) * phi * (RHS[lim(i - 1, N)][j] - 2.0f * RHS[i][j] + RHS[lim(i + 1, N)][j]);
                    d[j] = RHS[i][j] + diff_term + 0.5f * partRHS[i][j];
                }

                tridiag(la, lb, lc, c_prime, d_prime, N, d, result);
                for (int j = 0; j < N; j++)
                {
                    V[i][j] = result[j];
                }
            }

            // Update time step counter
            timeStepCounter++;
        }
    }
#endif // LINMONO || MONODOMAIN

#else
    if (strcmp(method, "SSI-ADI") == 0)
    {
        while (timeStepCounter < M)
        {
            // Get time step
            actualTime = time[timeStepCounter];

            // ================================================!
            //  Calcula Approx.                                !
            // ================================================!
            real x;
            for (int i = 0; i < N; i++)
            {
                x = i * delta_x;
                real diff_term = (sigma / (chi * Cm)) * 0.5f * phi * (V[lim(i - 1, N)] - 2.0f * V[i] + V[lim(i + 1, N)]);

#ifdef AFHN
                real RHS_V_term = RHS_V(V[i], W[i]);
#endif // AFHN

                Vtilde[i] = V[i] + diff_term + (0.5f * delta_t * (- RHS_V_term));

// Preparing part of the RHS of the following linear systems
#ifdef AFHN
                // Calculate W approximation
                real RHS_Vtilde_term = RHS_V(Vtilde[i], W[i]);
#endif // AFHN

                partRHS[i] = delta_t * (- RHS_Vtilde_term);

#ifdef AFHN
                // Update Wn+1
                real RHS_W_term = RHS_W(V[i], W[i]);
                real Wtilde = W[i] + (0.5f * delta_t * RHS_W_term);
                W[i] = W[i] + delta_t * RHS_W(Vtilde[i], Wtilde);
#endif // AFHN
            }

            // ================================================!
            //  Calcula V em n + 1 -> Resultado vai para V     !
            // ================================================!
            for (int i = 0; i < N; i++)
            {
                d[i] = V[i] + partRHS[i];
            }

            tridiag(la, lb, lc, c_prime, d_prime, N, d, result);
            for (int i = 0; i < N; i++)
            {
                V[i] = result[i];
            }

            // Update time step counter
            timeStepCounter++;
        }
    }
#endif // not CABLEEQ

    real finishTime = omp_get_wtime();
    real elapsedTime = finishTime - startTime;

// Calculate error
#ifdef CONVERGENCE_ANALYSIS
    real norm2error = calculateNorm2Error(V, exact, N, totalTime, delta_x);
#endif // CONVERGENCE_ANALYSIS

    // Write infos to file
    char infosFilePath[MAX_STRING_SIZE];
    snprintf(infosFilePath, MAX_STRING_SIZE * sizeof(char), "%s/infos/infos_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);
    FILE *fpInfos = fopen(infosFilePath, "w");
    printf("Infos saved to %s\n", infosFilePath);
    fprintf(fpInfos, "Domain Length = %d, Time = %f\n", L, totalTime);
    #ifndef CABLEEQ
    fprintf(fpInfos, "delta_x = %lf, Space steps N = %d, N*N = %d\n", delta_x, N, N * N);
    #else
    fprintf(fpInfos, "delta_x = %lf, Space steps N = %d\n", delta_x, N);
    #endif // not CABLEEQ
    fprintf(fpInfos, "delta_t = %lf, Time steps = %d\n", delta_t, M);
    fprintf(fpInfos, "Method %s\n", method);
    fprintf(fpInfos, "\nTotal execution time = %lf\n", elapsedTime);
#ifdef CONVERGENCE_ANALYSIS
    fprintf(fpInfos, "\nNorm-2 Error = %lf\n", norm2error);
#endif // CONVERGENCE_ANALYSIS
    fclose(fpInfos);

    // Save last frame
    char lastFrameFilePath[MAX_STRING_SIZE];
    snprintf(lastFrameFilePath, MAX_STRING_SIZE * sizeof(char), "%s/lastframe/last_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);
    FILE *fpLast = fopen(lastFrameFilePath, "w");
    printf("Last frame saved to %s\n", lastFrameFilePath);
#ifdef CONVERGENCE_ANALYSIS
    char exactFilePath[MAX_STRING_SIZE];
    snprintf(exactFilePath, MAX_STRING_SIZE * sizeof(char), "%s/exact/exact_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);
    FILE *fpExact = fopen(exactFilePath, "w");
    printf("Exact solution saved to %s\n", exactFilePath);
    char errorsFilePath[MAX_STRING_SIZE];
    snprintf(errorsFilePath, MAX_STRING_SIZE * sizeof(char), "%s/errors/errors_%.5f_%.5f.txt", pathToSaveData, delta_t, delta_x);
    FILE *fpErrors = fopen(errorsFilePath, "w");
    printf("Errors saved to %s\n", errorsFilePath);
#endif // CONVERGENCE_ANALYSIS
#ifndef CABLEEQ
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            fprintf(fpLast, "%e ", V[i][j]);
#ifdef CONVERGENCE_ANALYSIS
            fprintf(fpExact, "%e ", exact[i][j]);
            fprintf(fpErrors, "%e ", abs(V[i][j] - exact[i][j]));
#endif // CONVERGENCE_ANALYSIS
        }
        fprintf(fpLast, "\n");
#ifdef CONVERGENCE_ANALYSIS
        fprintf(fpExact, "\n");
        fprintf(fpErrors, "\n");
#endif // CONVERGENCE_ANALYSIS
    }
#else
    for (int i = 0; i < N; i++)
    {
        fprintf(fpLast, "%e ", V[i]);
    }
#endif // not CABLEEQ
    fclose(fpLast);
#ifdef CONVERGENCE_ANALYSIS
    fclose(fpExact);
    fclose(fpErrors);
#endif // CONVERGENCE_ANALYSIS

    // Free memory
    free(time);

#ifndef CABLEEQ
    for (int i = 0; i < N; i++)
    {
        free(V[i]);
        free(Vtilde[i]);
        free(RHS[i]);
        free(partRHS[i]);
        free(exact[i]);
#ifdef MONODOMAIN
#ifdef AFHN
        free(W[i]);
#endif // AFHN
#endif // MONODOMAIN
    }
#endif // not CABLEEQ
    free(V);
    free(Vtilde);
    free(RHS);
    free(partRHS);
    free(exact);
    free(c_prime);
    free(d_prime);
    free(d);
    free(result);
    free(la);
    free(lb);
    free(lc);
    free(pathToSaveData);
#if defined(MONODOMAIN) || defined(CABLEEQ)
#ifdef AFHN
    free(W);
#endif // AFHN
#endif // MONODOMAIN || CABLEEQ

    return;
}

#endif // CPU_METHODS_H