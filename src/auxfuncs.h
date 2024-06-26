#ifndef AUXFUNCS_H
#define AUXFUNCS_H

#include "include.h"
#include "functions.h"

// Populate diagonals for Thomas algorithm
void populateDiagonalThomasAlgorithm(real* la, real* lb, real* lc, int N, real phi)
{
    for (int i = 0; i < N; i++)
    {
        la[i] = -phi;
        lb[i] = 1.0 + 2.0*phi;
        lc[i] = -phi;
    }
    lc[0] = lc[0] + la[0];
    la[N-1] = la[N-1] + lc[N-1];
    la[0] = 0.0;
    lc[N-1] = 0.0; 
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
void initializeVariableWithExactSolution(real** Var, int N, real delta_x)
{
    real x, y;
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            x = i * delta_x;
            y = j * delta_x;
            Var[i][j] = exactSolution(0.0, x, y);
        }
    }
}

void initializeVariableWithValue(real** Var, int N, real delta_x, real value)
{
    real x, y;
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            x = i * delta_x;
            y = j * delta_x;
            Var[i][j] = value;
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
                // diffusion = (V[i+1][j] - V[i][j]);
            else if (i == N-1)
                diffusion = (2.0*V[i-1][j] - 2.0*V[i][j]);
                // diffusion = (V[i-1][j] - V[i][j]);

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
                // diffusion = (V[i][j+1] - V[i][j]);
            else if (j == N-1)
                diffusion = (2.0*V[i][j-1] - 2.0*V[i][j]);
                // diffusion = (V[i][j-1] - V[i][j]);

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

int lim(int num, int N)
{
    if (num == -1)
        return 1;
    else if (num == N)
    {
        return N-2;
    }
    return num;
}

real calculateNorm2Error(real** V, real** exact, int N, real T, real delta_x)
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
            exact[i][j] = solution;
            sum += ((V[i][j] - solution) * (V[i][j] - solution));
        }
    }
    return delta_x * sqrt(sum);
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

#endif // SERIAL

#endif // AUXFUNCS_H