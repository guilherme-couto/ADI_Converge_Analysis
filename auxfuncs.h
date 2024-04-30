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

void initializeTimeArray(real* timeArray, int M, real dt)
{
    for (int i = 0; i < M; ++i)
    {
        timeArray[i] = i * dt;
    }
}

#endif // AUXFUNCS_H