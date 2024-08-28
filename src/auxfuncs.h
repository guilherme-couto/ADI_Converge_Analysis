#ifndef AUXFUNCS_H
#define AUXFUNCS_H

#ifdef GPU
#include "gpu_functions.h"
#else
#include "cpu_functions.h"
#endif

// Populate diagonals for Thomas algorithm
// la -> Subdiagonal
// lb -> Diagonal
// lc -> Superdiagonal
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

void createDirectories(char* method, real theta, char* pathToSaveData)
{   
    // Build the path
    char path[MAX_STRING_SIZE];
    snprintf(path, MAX_STRING_SIZE*sizeof(char), "./simulation_files/%s/AFHN/%s", REAL_TYPE, method);

    // Add theta to the path
    if (strcmp(method, "theta-ADI") == 0) {
        char thetaPath[MAX_STRING_SIZE];
        snprintf(thetaPath, MAX_STRING_SIZE*sizeof(char), "%.2lf", theta);
        strcat(path, "/");
        strcat(path, thetaPath);
    }

    // Make directories
    char command[MAX_STRING_SIZE];
    snprintf(command, MAX_STRING_SIZE*sizeof(char), "mkdir -p %s", path);
    system(command);

    // Update pathToSaveData
    strcpy(pathToSaveData, path);
}

void initializeTimeArray(real *timeArray, int M, real dt)
{
    for (int i = 0; i < M; ++i)
    {
        timeArray[i] = i * dt;
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

#ifdef MONOAFHN
void populateStimuli(Stimulus *stimuli, real delta_x)
{
    for (int i = 0; i < numberOfStimuli; i++)
    {
        stimuli[i].strength = stimuliStrength;
        stimuli[i].begin = stimuliBegin[i];
        stimuli[i].duration = stimuliDuration;
        
        // Discretized limits of stimulation areas
        stimuli[i].xMaxDisc = round(stimulixMax[i] / delta_x);
        stimuli[i].xMinDisc = round(stimulixMin[i] / delta_x);
        stimuli[i].yMaxDisc = round(stimuliyMax[i] / delta_x);
        stimuli[i].yMinDisc = round(stimuliyMin[i] / delta_x);
    }
}
#endif // MONOAFHN

#ifdef SERIAL
void initialize2DVariableWithExactSolution(real** Var, int N, real delta_x)
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

void initialize2DVariableWithValue(real** Var, int N, real value)
{
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            Var[i][j] = value;
        }
    }
}

void initialize2DVariableFromFile(real** Var, int N, char* filename, real delta_x)
{
    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        printf("Error opening file %s\n", filename);
        exit(1);
    }

    printf("Reading file %s to initialize variable\n", filename);
    int baseN = round(L / 0.0005) + 1;
    int rate = round(delta_x / 0.0005);
    real value;
    for (int i = 0; i < baseN; ++i)
    {
        for (int j = 0; j < baseN; ++j)
        {
            // Read value from file to variable
            // If i and j are multiples of rate, read value to Var
            #ifndef USE_DOUBLE
            fscanf(file, "%e", &value);
            #else
            fscanf(file, "%le", &value);
            #endif
            if (i % rate == 0 && j % rate == 0)
            {
                Var[(i/rate)][(j/rate)] = value;
            }
        }
    }
    fclose(file);
    printf("Variable initialized with values from file %s\n", filename);
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

void saveFrame(FILE *file, real actualTime, real** V, int N)
{
    fprintf(file, "%lf\n", actualTime);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            fprintf(file, "%e ", V[i][j]);
        }
        fprintf(file, "\n");
    }
}
#endif // SERIAL

#ifdef GPU
void initialize2DVariableWithValue(real* Var, int N, real value)
{
    int index;
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            index = i*N + j;
            Var[index] = value;
        }
    }
}

void initialize2DVariableFromFile(real* Var, int N, char* filename, real delta_x)
{
    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        printf("Error opening file %s\n", filename);
        exit(1);
    }

    printf("Reading file %s to initialize variable\n", filename);
    int index;
    int baseN = round(L / 0.0005) + 1;
    int rate = round(delta_x / 0.0005);
    real value;
    for (int i = 0; i < baseN; ++i)
    {
        for (int j = 0; j < baseN; ++j)
        {
            // Read value from file to variable
            // If i and j are multiples of rate, read value to Var
            #ifndef USE_DOUBLE
            fscanf(file, "%e", &value);
            #else
            fscanf(file, "%le", &value);
            #endif
            if (i % rate == 0 && j % rate == 0)
            {
                index = (i/rate)*N + (j/rate);
                Var[index] = value;
            }
        }
    }
    fclose(file);
    printf("Variable initialized with values from file %s\n", filename);
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

void saveFrame(FILE *file, real actualTime, real* V, int N)
{
    fprintf(file, "%lf\n", actualTime);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            int index = i*N + j;
            fprintf(file, "%e ", V[index]);
        }
        fprintf(file, "\n");
    }
}

#ifdef CONVERGENCE_ANALYSIS
void initialize2DVariableWithExactSolution(real* Var, int N, real delta_x)
{
    real x, y;
    int index;
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            x = i * delta_x;
            y = j * delta_x;
            index = i*N + j;
            Var[index] = exactSolution(0.0, x, y);
        }
    }
}

real calculateNorm2Error(real* V, real** exact, int N, real T, real delta_x)
{
    real x, y;
    int index;
    real solution;
    real sum = 0.0;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            x = j * delta_x;
            y = i * delta_x;
            index = i*N + j;
            solution = exactSolution(T, x, y);
            exact[i][j] = solution;
            sum += ((V[index] - solution) * (V[index] - solution));
        }
    }
    return delta_x * sqrt(sum);
}
#endif // CONVERGENCE_ANALYSIS
#endif // GPU

#endif // AUXFUNCS_H