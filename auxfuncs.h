#ifndef AUXFUNCS_H
#define AUXFUNCS_H

#include "include.h"

// Populate diagonals for Thomas algorithm
void populateDiagonalThomasAlgorithm(real* la, real* lb, real* lc, int N, real phi)
{
    // First row
    la[0] = 0.0;
    real b = 1 + 2 * phi;     // diagonal (1st and last row)
    real c = - 2 * phi;       // superdiagonal
    lb[0] = b;
    lc[0] = c;

    // Middle rows
    real a = -phi;            // subdiagonal
    c = -phi;                   // superdiagonal
    for (int i = 1; i < N - 1; ++i)
    {
        la[i] = a;
        lb[i] = b;
        lc[i] = c;
    }

    // Last row
    a = - 2 * phi;
    la[N - 1] = a;
    lb[N - 1] = b;
    lc[N - 1] = 0.0;
}

void createDirectoriesAndFiles(char* method, real theta, char* pathToSaveData, char* aux)
{
    char command[MAX_STRING_SIZE];
    sprintf(command, "%s", "mkdir -p");
    char path[MAX_STRING_SIZE] = "./simulation-files/";
    strcat(path, REAL_TYPE);
    sprintf(aux, "%s %s", command, path);
    system(aux);
    sprintf(pathToSaveData, "%s/%s", path, "AFHN");
    sprintf(aux, "%s %s", command, pathToSaveData);
    system(aux);
    sprintf(pathToSaveData, "%s/%s", pathToSaveData, method);
    sprintf(aux, "%s %s", command, pathToSaveData);
    system(aux);
    if (strcmp(method, "theta-ADI") == 0)
    {
        sprintf(pathToSaveData, "%s/%.2lf", pathToSaveData, theta);
        sprintf(aux, "%s %s", command, pathToSaveData);
        system(aux);
    }
}

void initializeTimeArray(real* timeArray, int M, real dt)
{
    for (int i = 0; i < M; ++i)
    {
        timeArray[i] = i * dt;
    }
}

#ifdef PARALLEL
void initializeStateVariable(real* V, int N)
{
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            V[i * N + j] = V_init;
        }
    }
}

void thomasFactorConstantBatch(real* la, real* lb, real* lc, int n) {

	int rowCurrent;
	int rowPrevious;

	rowCurrent = 0;

	// First row
	lb[rowCurrent] = lb[rowCurrent];
	lc[rowCurrent] = lc[rowCurrent] / lb[rowCurrent];

	for (int i = 1; i < n - 1; ++i)	{
		rowPrevious = rowCurrent;
		rowCurrent  += 1;

		la[rowCurrent] = la[rowCurrent];
		lb[rowCurrent] = lb[rowCurrent] - la[rowCurrent]*lc[rowPrevious];
		lc[rowCurrent] = lc[rowCurrent] / lb[rowCurrent];
	}

	rowPrevious = rowCurrent;
	rowCurrent += 1;

	// Last row
	la[rowCurrent] = la[rowCurrent];
	lb[rowCurrent] = lb[rowCurrent] - la[rowCurrent]*lc[rowPrevious];
}
#endif // PARALLEL

#ifdef SERIAL
void initializeStateVariable(real** V, int N)
{
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            V[i][j] = V_init;
        }
    }
}

real forcingTerm(real x, real y, real t)
{
    return cos(_pi*x/L) * cos(_pi*y/L) * (chi*Cm*exp(-t) + ((2.0*_pi*_pi*sigma)/(L*L))*(1.0-exp(-t)) + (chi*G)*(1.0-exp(-t))); 
}

void calculateVApprox(real** V, real** Rv, int N, real delta_x, real delta_t, real actual_t)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            real actualV = V[i][j];

            // Diffusion component
            real diffusion = 0.0;

            // Diffusion in x
            if (j > 0 && j < N-1)
                diffusion = (V[i][j+1] - 2.0*V[i][j] + V[i][j-1]) / (delta_x*delta_x);
            else if (j == 0)
                diffusion = (2.0*V[i][j+1] - 2.0*V[i][j]) / (delta_x*delta_x);
            else if (j == N-1)
                diffusion = (2.0*V[i][j-1] - 2.0*V[i][j]) / (delta_x*delta_x);

            // Diffusion in y
            if (i > 0 && i < N-1)
                diffusion += (V[i+1][j] - 2.0*V[i][j] + V[i-1][j]) / (delta_x*delta_x);
            else if (i == 0)
                diffusion += (2.0*V[i+1][j] - 2.0*V[i][j]) / (delta_x*delta_x);
            else if (i == N-1)
                diffusion += (2.0*V[i-1][j] - 2.0*V[i][j]) / (delta_x*delta_x);

            // Forcing term
            real x = j * delta_x;
            real y = i * delta_x;
            real forcing = forcingTerm(x, y, actual_t);

            // Calculate the RHS with actual values
            real actualRHS = (forcing/(chi*Cm)) - (G*actualV/Cm);

            // Calculate an approx for V
            real Vtilde;
            Vtilde = actualV + (delta_t/2.0) * (((sigma/(chi*Cm)) * diffusion) + actualRHS);

            // Recalculate the forcing term at time t+(dt/2) and the RHS with the new values
            real new_t = actual_t + (0.5 * delta_t);
            forcing = forcingTerm(x, y, new_t);
            real tildeRHS = (forcing/(chi*Cm)) - (G*Vtilde/Cm);

            // Update reaction term
            Rv[i][j] = tildeRHS;
        }
    }
}

void prepareRHS(real** V, real** RHS, real** Rv, int N, real phi, real delta_t)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            // Explicit diffusion in y
            real diffusion = 0.0;
            if (i > 0 && i < N-1)
                diffusion = (V[i+1][j] - 2.0*V[i][j] + V[i-1][j]);
            else if (i == 0)
                diffusion = (2.0*V[i+1][j] - 2.0*V[i][j]);
            else if (i == N-1)
                diffusion = (2.0*V[i-1][j] - 2.0*V[i][j]);

            RHS[i][j] = V[i][j] + (phi * diffusion) + (0.5*delta_t*Rv[i][j]);
        }
    }
}
#endif // SERIAL

#endif // AUXFUNCS_H