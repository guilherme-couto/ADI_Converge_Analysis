#ifndef AUXFUNCS_H
#define AUXFUNCS_H

#include "include.h"

// Populate diagonals for Thomas algorithm
void populateDiagonalThomasAlgorithm(real* la, real* lb, real* lc, int N, real phi)
{
    // First row
    la[0] = 0.0;
    real b = 1.0 + 2.0 * phi;     // diagonal (1st and last row)
    real c = - 2.0 * phi;       // superdiagonal
    lb[0] = b;
    lc[0] = c;

    // Middle rows
    real a = -phi;            // subdiagonal
    c = -phi;                 // superdiagonal
    for (int i = 1; i < N - 1; ++i)
    {
        la[i] = a;
        lb[i] = b;
        lc[i] = c;
    }

    // Last row
    a = - 2.0 * phi;
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
void initializeStateVariable(real** V, int N, real delta_x)
{
    real x, y;
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            x = j * delta_x;
            y = i * delta_x;
            V[i][j] = exactSolution(0.0, x, y);
        }
    }
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
            real forcing = forcingTerm(x, y, actual_t+(delta_t*0.5));

            // Aux variables
            real actualRHS, Vtilde, tildeRHS;

            #ifdef LINMONO
            // Calculate the RHS with actual values
            actualRHS = (forcing/(chi*Cm)) - (G*actualV/Cm);

            // Calculate an approx for V
            Vtilde = actualV + (0.5 * delta_t) * (((sigma/(chi*Cm)) * diffusion) + actualRHS);

            // Recalculate the forcing term at time t+(dt/2) and the RHS with the new values
            tildeRHS = -(G*Vtilde/Cm);
            #endif // LINMONO

            // Update reaction term
            Rv[i][j] = tildeRHS;
        }
    }
}

void prepareRHS_explicit_y(real** V, real** RHS, real** Rv, int N, real phi, real delta_t, real time, real delta_x)
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

            real x = j * delta_x;
            real y = i * delta_x;

            #ifdef LINMONO
            RHS[i][j] = V[i][j] + (phi * diffusion) + (0.5*delta_t*Rv[i][j]) + (0.5*delta_t*forcingTerm(x, y, time)/(chi*Cm));
            #endif // LINMONO
            #ifdef DIFF
            RHS[i][j] = V[i][j] + (phi * diffusion) + (0.5*delta_t*forcingTerm(x, y, time));
            #endif // DIFF
        }
    }
}

void prepareRHS_explicit_x(real** V, real** RHS, real** Rv, int N, real phi, real delta_t, real time, real delta_x)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            // Explicit diffusion in x
            real diffusion = 0.0;
            if (j > 0 && j < N-1)
                diffusion = (V[i][j+1] - 2.0*V[i][j] + V[i][j-1]);
            else if (j == 0)
                diffusion = (2.0*V[i][j+1] - 2.0*V[i][j]);
            else if (j == N-1)
                diffusion = (2.0*V[i][j-1] - 2.0*V[i][j]);

            real x = j * delta_x;
            real y = i * delta_x;
            
            #ifdef LINMONO
            RHS[i][j] = V[i][j] + (phi * diffusion) + (0.5*delta_t*Rv[i][j]) + (0.5*delta_t*forcingTerm(x, y, time)/(chi*Cm));
            #endif // LINMONO
            #ifdef DIFF
            RHS[i][j] = V[i][j] + (phi * diffusion) + (0.5*delta_t*forcingTerm(x, y, time));
            #endif // DIFF
        }
    }
}

real calculateNorm2Error(real**V, int N, real T, real delta_x)
{
    real x, y;
    real solution;
    real sum = 0.0;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            x = j * delta_x;
            y = i * delta_x;
            solution = exactSolution(T, x, y);
            sum += (V[i][j] - solution) * (V[i][j] - solution);
        }
    }
    return sqrt(delta_x*delta_x*sum);
}

void copyMatrices(real** in, real** out, int N)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            out[i][j] = in[i][j];
        }
    }
}

void copyColumnToVector(real** V, real* d, int N, int column_id)
{
    for (int i = 0; i < N; i++)
    {
        d[i] = V[i][column_id];
    }
}

void copyVectorToColumn(real** V, real* d, int N, int column_id)
{
    for (int i = 0; i < N; i++)
    {
        V[i][column_id] = d[i];
    }
}

void thomasAlgorithm(real* la, real* lb, real* lc, real* c_prime, real* d_prime, int N, real* d)
{
    // la -> subdiagonal
    // lb -> diagonal
    // lc -> superdiagonal
    // d -> RHS initially and, in the end, will store the result

    // 1st: Forward sweep
    c_prime[0] = lc[0] / lb[0];
    d_prime[0] = d[0] / lb[0];
    for (int i = 1; i < N; i++)
    {
        if (i < N-1)
            c_prime[i] = lc[i] / (lb[i] - (la[i]*c_prime[i-1]));
        
        d_prime[i] = (d[i] - (la[i]*d_prime[i-1])) / (lb[i] - (la[i]*c_prime[i-1]));
    }

    // 2nd: Back substitution
    d[N-1] = d_prime[N-1];
    for (int i = N-2; i >= 0; i--)
    {
        d[i] = d_prime[i] - (c_prime[i]*d[i+1]);
    }

    // Vector d now has the result
}

// Tridiag Ricardo
void tridiag(real* la, real* lb, real* lc, real* c_prime, real* d_prime, int N, real* d, real* result)
{
    c_prime[0] = lc[0]/lb[0];
    for (int i = 1; i < N-1; i++)
    {
        c_prime[i] = lc[i] / (lb[i] - c_prime[i-1]*la[i]);
    }
    d_prime[0] = d[0]/lb[0];
    for (int i = 1; i < N; i++)
    {
        d_prime[i] = (d[i]-d_prime[i-1]*la[i]) / (lb[i]-c_prime[i-1]*la[i]);
    }
    result[N-1] = d_prime[N-1];
    for (int i = N-2; i >= 0; i--)
    {
        result[i] = d_prime[i] - c_prime[i]*result[i+1];
    }
}

void solveExplicitly(real** V, real** RHS, int N, real delta_x, real delta_t, real time)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            real actualV = V[i][j];
            // if (time < delta_t && i == 0 && j == 0)
            // {
            //     printf("%f/%lf ", time, actualV);
            // }
                

            // Diffusion component
            real diffusion_x = 0.0;
            real diffusion_y = 0.0;

            // Diffusion in x
            if (j > 0 && j < N-1)
                diffusion_x = (V[i][j+1] - 2.0*actualV + V[i][j-1]) / (delta_x*delta_x);
            else if (j == 0)
                diffusion_x = (2.0*V[i][j+1] - 2.0*actualV) / (delta_x*delta_x);
            else if (j == N-1)
                diffusion_x = (2.0*V[i][j-1] - 2.0*actualV) / (delta_x*delta_x);
            
            // Diffusion in y
            if (i > 0 && i < N-1)
                diffusion_y = (V[i+1][j] - 2.0*actualV + V[i-1][j]) / (delta_x*delta_x);
            else if (i == 0)
                diffusion_y = (2.0*V[i+1][j] - 2.0*actualV) / (delta_x*delta_x);
            else if (i == N-1)
                diffusion_y = (2.0*V[i-1][j] - 2.0*actualV) / (delta_x*delta_x);

            // Forcing term
            real x = j * delta_x;
            real y = i * delta_x;
            real forcing = forcingTerm(x, y, time);

            // Calculate new value of V
            #ifdef LINMONO
            RHS[i][j] = actualV + delta_t * (((sigma/(chi*Cm)) * (diffusion_x + diffusion_y)) + (forcing/(chi*Cm)) - (G*actualV/Cm));
            #endif // LINMONO
            #ifdef DIFF
            RHS[i][j] = actualV + delta_t * ((sigma * (diffusion_x + diffusion_y)) + forcing);
            #endif // DIFF
            if (i == 0 && j == 0)
            {
                printf("time %lf - RHS[0][0] = %e\n", time, RHS[0][0]);
            }
            if (time == 0 && i == 0 && j == 0)
            {
                printf("diff_x = %e\n", diffusion_x);
            }
        }
    }

    // Copy from RHS to V
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            V[i][j] = RHS[i][j];
}
#endif // SERIAL

#endif // AUXFUNCS_H