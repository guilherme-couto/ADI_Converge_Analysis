#ifndef AUXFUNCS_H
#define AUXFUNCS_H

#ifdef GPU
#include "gpu_functions.h"
#else
#include "cpu_functions.h"
#endif

#ifdef SERIAL
void initialize2DVariableWithExactSolution(real **Var, int Nx, int Ny, real delta_x, real delta_y)
{
    real x, y;
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            y = i * delta_y;
            x = j * delta_x;
            Var[i][j] = exactSolution(0.0f, x, y);
        }
    }
}

void initialize1DVariableWithValue(real *Var, int N, real value)
{
    for (int i = 0; i < N; i++)
    {
        Var[i] = value;
    }   
}

void initialize2DVariableWithValue(real **Var, int Nx, int Ny, real value)
{
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            Var[i][j] = value;
        }
    }
}

void initialize1DVariableFromFile(real *Var, int N, char *filename, real delta_x, char *varName, real reference_dx)
{
    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        ERRORMSG("Error opening file %s\n", filename);
        exit(1);
    }

    int baseN = round(Lx / reference_dx) + 1;
    int rate = round(delta_x / reference_dx);

    int sizeFile = 0;
    int sizeVar = 0;
    real value;

    INFOMSG("Reading file %s to initialize variable with a rate of %d\n", filename, rate);

    for (int i = 0; i < baseN; i++)
    {
        // Read value from file to variable
        // If i is multiple of rate, read value to Var
        fscanf(file, FSCANF_REAL, &value);
        if (i % rate == 0)
        {
            Var[int(i / rate)] = value;
            if (isnan(value))
            {
                ERRORMSG("At var index [%d], file index %d, value is NaN\n", int(i / rate), i);
                exit(1);
            }
            sizeVar++;
        }
        sizeFile++;
        
    }
    fclose(file);

    SUCCESSMSG("Variable %s initialized with %d values from the %d values in file\n", varName, sizeVar, sizeFile);
}

void initialize2DVariableFromFile(real **Var, int Nx, int Ny, char *filename, real delta_x, real delta_y, char *varName, real reference_dx, real reference_dy)
{
    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        ERRORMSG("Error opening file %s\n", filename);
        exit(1);
    }

    int baseNx = round(Lx / reference_dx) + 1;
    int baseNy = round(Ly / reference_dy) + 1;
    int rate_x = round(delta_x / reference_dx);
    int rate_y = round(delta_y / reference_dy);

    int sizeFile = 0;
    int sizeVar = 0;
    real value;

    INFOMSG("Reading file %s to initialize variable with rate_x=%d and rate_y=%d\n", filename, rate_x, rate_y);
    for (int i = 0; i < baseNy; i++)
    {
        for (int j = 0; j < baseNx; j++)
        {
            // Read value from file to variable
            // If i and j are multiples of rate, read value to Var
            fscanf(file, FSCANF_REAL, &value);
            if (i % rate_y == 0 && j % rate_x == 0)
            {
                Var[int(i / rate_y)][int(j / rate_x)] = value;
                if (isnan(value))
                {
                    ERRORMSG("At var index [%d][%d], file index %d, value is NaN\n", int(i / rate_y), int(j / rate_x), i * baseNx + j);
                    exit(1);
                }
                sizeVar++;
            }
            sizeFile++;
        }
    }
    fclose(file);
    INFOMSG("Variable %s initialized with %d values from the %d values in file\n", varName, sizeVar, sizeFile);
}

void shift1DVariableToLeft(real *Var, int N, real length, real delta_x, real initValue, char *varName)
{
    real *temp = (real *)malloc(N * sizeof(real));
    for (int i = 0; i < N; i++)
    {
        temp[i] = Var[i];
    }
    
    int lengthIndex = round(length / delta_x) + 1;
    for (int i = 0; i < N - lengthIndex; i++)
    {
        Var[i] = temp[i + lengthIndex];
    }
    for (int i = N - lengthIndex; i < N; i++)
    {
        Var[i] = initValue;
    }
    free(temp);
    SUCCESSMSG("Variable %s shifted to the left by %.2f cm\n", varName, length);
}

void shift2DVariableToLeft(real **Var, int Nx, int Ny, real length, real delta_x, real delta_y, real initValue, char *varName)
{
    real *temp = (real *)malloc(Nx * sizeof(real));
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            temp[j] = Var[i][j];
        }

        int lengthIndex = round(length / delta_x) + 1;
        for (int j = 0; j < Nx - lengthIndex; j++)
        {
            Var[i][j] = temp[j + lengthIndex];
        }
        for (int j = Nx - lengthIndex; j < Nx; j++)
        {
            Var[i][j] = initValue;
        }
    }
    free(temp);
    SUCCESSMSG("Variable %s shifted to the left by %.2f cm\n", varName, length);
}

real calculateNorm2Error(real **V, real **exact, int Nx, int Ny, real totalTime, real delta_x, real delta_y)
{
    real x, y;
    real solution;
    real sum = 0.0;
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            x = j * delta_x;
            y = i * delta_y;
            solution = exactSolution(totalTime, x, y);
            exact[i][j] = solution;
            sum += ((V[i][j] - solution) * (V[i][j] - solution));
        }
    }
    return sqrt(sum / (Nx * Ny));
}

void copyMatrices(real **in, real **out, int Nx, int Ny)
{
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            out[i][j] = in[i][j];
        }
    }
}

void copyColumnToVector(real **V, real *d, int N, int column_id)
{
    for (int i = 0; i < N; i++)
    {
        d[i] = V[i][column_id];
    }
}

void copyVectorToColumn(real **V, real *d, int N, int column_id)
{
    for (int i = 0; i < N; i++)
    {
        V[i][column_id] = d[i];
    }
}

void thomasAlgorithm(real *la, real *lb, real *lc, real *c_prime, real *d_prime, int N, real *d)
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
        if (i < N - 1)
            c_prime[i] = lc[i] / (lb[i] - (la[i] * c_prime[i - 1]));

        d_prime[i] = (d[i] - (la[i] * d_prime[i - 1])) / (lb[i] - (la[i] * c_prime[i - 1]));
    }

    // 2nd: Back substitution
    d[N - 1] = d_prime[N - 1];
    for (int i = N - 2; i >= 0; i--)
    {
        d[i] = d_prime[i] - (c_prime[i] * d[i + 1]);
    }

    // Vector d now has the result
}

void tridiag(real *la, real *lb, real *lc, real *c_prime, real *d_prime, int N, real *d, real *result)
{
    c_prime[0] = lc[0] / lb[0];
    for (int i = 1; i < N - 1; i++)
    {
        c_prime[i] = lc[i] / (lb[i] - c_prime[i - 1] * la[i]);
    }
    d_prime[0] = d[0] / lb[0];
    for (int i = 1; i < N; i++)
    {
        d_prime[i] = (d[i] - d_prime[i - 1] * la[i]) / (lb[i] - c_prime[i - 1] * la[i]);
    }
    result[N - 1] = d_prime[N - 1];
    for (int i = N - 2; i >= 0; i--)
    {
        result[i] = d_prime[i] - c_prime[i] * result[i + 1];
    }
}

#ifndef CABLEEQ
void saveFrame(FILE *file, real actualTime, real **V, int Nx, int Ny)
{
    fprintf(file, "%lf\n", actualTime);
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            fprintf(file, "%e ", V[i][j]);
        }
        fprintf(file, "\n");
    }
}
#else
void saveFrame(FILE *file, real actualTime, real *V, int N)
{
    fprintf(file, "%lf\n", actualTime);
    for (int i = 0; i < N; i++)
    {
        fprintf(file, "%e ", V[i]);
    }
    fprintf(file, "\n");
}
#endif // not CABLEEQ
#endif // SERIAL

#ifdef GPU
void initialize2DVariableWithValue(real *Var, int Nx, int Ny, real value)
{
    int index;
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            index = i * Nx + j;
            Var[index] = value;
        }
    }
}

void initialize2DVariableFromFile(real *Var, int Nx, int Ny, char *filename, real delta_x, real delta_y, char *varName, real reference_dx, real reference_dy)
{
    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        ERRORMSG("Error opening file %s\n", filename);
        exit(1);
    }

    int baseNx = round(Lx / reference_dx) + 1;
    int baseNy = round(Ly / reference_dy) + 1;
    int rate_x = round(delta_x / reference_dx);
    int rate_y = round(delta_y / reference_dy);

    int sizeFile = 0;
    int sizeVar = 0;
    real value;

    INFOMSG("Reading file %s to initialize variable with rate_x=%d and rate_y=%d\n", filename, rate_x, rate_y);
    for (int i = 0; i < baseNy; i++)
    {
        for (int j = 0; j < baseNx; j++)
        {
            // Read value from file to variable
            // If i and j are multiples of rate, read value to Var
            fscanf(file, FSCANF_REAL, &value);
            if (i % rate_y == 0 && j % rate_x == 0)
            {
                Var[sizeVar] = value;
                if (isnan(value))
                {
                    ERRORMSG("At var index %d, file index %d, value is NaN\n", sizeVar, i * baseNx + j);
                    exit(1);
                }
                sizeVar++;
            }
            sizeFile++;
        }
    }
    fclose(file);
    INFOMSG("Variable %s initialized with %d values from the %d values in file\n", varName, sizeVar, sizeFile);
}

void shift2DVariableToLeft(real *Var, int N, real length, real delta_x, real initValue, char *varName)
{
    real *temp = (real *)malloc(N * sizeof(real));
    int index;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            index = i * N + j;
            temp[j] = Var[index];
        }

        int lengthIndex = round(length / delta_x) + 1;
        for (int j = 0; j < N - lengthIndex; j++)
        {
            index = i * N + j;
            Var[index] = temp[j + lengthIndex];
        }
        for (int j = N - lengthIndex; j < N; j++)
        {
            index = i * N + j;
            Var[index] = initValue;
        }
    }
    free(temp);
    SUCCESSMSG("Variable %s shifted to the left by %.2f cm\n", varName, length);
}

void thomasFactorConstantBatch(real *la, real *lb, real *lc, int n)
{

    int rowCurrent;
    int rowPrevious;

    rowCurrent = 0;

    // First row
    lb[rowCurrent] = lb[rowCurrent];
    lc[rowCurrent] = lc[rowCurrent] / lb[rowCurrent];

    for (int i = 1; i < n - 1; i++)
    {
        rowPrevious = rowCurrent;
        rowCurrent += 1;

        la[rowCurrent] = la[rowCurrent];
        lb[rowCurrent] = lb[rowCurrent] - la[rowCurrent] * lc[rowPrevious];
        lc[rowCurrent] = lc[rowCurrent] / lb[rowCurrent];
    }

    rowPrevious = rowCurrent;
    rowCurrent += 1;

    // Last row
    la[rowCurrent] = la[rowCurrent];
    lb[rowCurrent] = lb[rowCurrent] - la[rowCurrent] * lc[rowPrevious];
}

void saveFrame(FILE *file, real actualTime, real *V, int N)
{
    fprintf(file, "%lf\n", actualTime);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            int index = i * N + j;
            fprintf(file, "%e ", V[index]);
        }
        fprintf(file, "\n");
    }
}

#ifdef CONVERGENCE_ANALYSIS_FORCING_TERM
void initialize2DVariableWithExactSolution(real *Var, int Nx, int Ny, real delta_x, real delta_y)
{
    real x, y;
    int index;
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            x = j * delta_x;
            y = i * delta_y;
            index = i * Nx + j;
            Var[index] = exactSolution(0.0f, x, y);
        }
    }
}

real calculateNorm2Error(real *V, real **exact, int N, real totalTime, real delta_x)
{
    real x, y;
    int index;
    real solution;
    real sum = 0.0f;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            x = j * delta_x;
            y = i * delta_x;
            index = i * N + j;
            solution = exactSolution(totalTime, x, y);
            exact[i][j] = solution;
            sum += ((V[index] - solution) * (V[index] - solution));
        }
    }
    return delta_x * sqrt(sum);
}
#endif // CONVERGENCE_ANALYSIS_FORCING_TERM
#endif // GPU

// Populate diagonals for Thomas algorithm
void populateDiagonalThomasAlgorithm(real *la, real *lb, real *lc, int N, real phi)
{
    // la -> Subdiagonal
    // lb -> Diagonal
    // lc -> Superdiagonal

    for (int i = 0; i < N; i++)
    {
        la[i] = -phi;
        lb[i] = 1.0f + 2.0f * phi;
        lc[i] = -phi;
    }
    lc[0] = lc[0] + la[0];
    la[N - 1] = la[N - 1] + lc[N - 1];
    la[0] = 0.0f;
    lc[N - 1] = 0.0f;
}

void createDirectories(real delta_t, real delta_x, real delta_y, char *pathToSaveData)
{
    // Build the path
    char path[MAX_STRING_SIZE];
    #ifndef CABLEEQ
    snprintf(path, MAX_STRING_SIZE * sizeof(char), "./simulation_files/dt_%.5g_dx_%.5g_dy_%.5g/%s/%s/%s/%s/%s", delta_t, delta_x, delta_y, EXECUTION_TYPE, REAL_TYPE, PROBLEM, CELL_MODEL, METHOD);
    #else
    snprintf(path, MAX_STRING_SIZE * sizeof(char), "./simulation_files/dt_%.5g_dx_%.5g/%s/%s/%s/%s/%s", delta_t, delta_x, EXECUTION_TYPE, REAL_TYPE, PROBLEM, CELL_MODEL, METHOD);
    #endif // not CABLEEQ

#ifdef THETA
    // Add theta to the path
    char thetaPath[MAX_STRING_SIZE];
    snprintf(thetaPath, MAX_STRING_SIZE * sizeof(char), "%.2lf", THETA);
    strcat(path, "/");
    strcat(path, thetaPath);
#endif // THETA

    // Update pathToSaveData
    strcpy(pathToSaveData, path);
}

void initializeTimeArray(real *timeArray, int M, real dt)
{
    for (int i = 0; i < M; i++)
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
        return N - 2;
    }
    return num;
}

#if defined(MONODOMAIN) || defined(CABLEEQ)
#ifndef CONVERGENCE_ANALYSIS_FORCING_TERM
void populateStimuli(Stimulus *stimuli, real delta_x, real delta_y)
{
    for (int i = 0; i < numberOfStimuli; i++)
    {
        stimuli[i].strength = stimuliStrength;
        stimuli[i].begin = stimuliBegin[i];
        stimuli[i].duration = stimuliDuration;

        // Discretized limits of stimulation areas
        stimuli[i].xMaxDisc = round(stimulixMax[i] / delta_x);
        stimuli[i].xMinDisc = round(stimulixMin[i] / delta_x);
        stimuli[i].yMaxDisc = round(stimuliyMax[i] / delta_y);
        stimuli[i].yMinDisc = round(stimuliyMin[i] / delta_y);

        #if defined(INIT_WITH_SPIRAL) || defined(RESTORE_STATE_AND_SHIFT)
        stimuli[i].strength = 0.0f;
        #endif // INIT_WITH_SPIRAL || RESTORE_STATE_AND_SHIFT

        #ifdef CABLEEQ
        // Only one stimulus for CABLEEQ
        if (i > 0)
            stimuli[i].strength = 0.0f;
        #endif // CABLEEQ
    }
}
#endif // not CONVERGENCE_ANALYSIS_FORCING_TERM
#endif // MONODOMAIN || CABLEEQ

#endif // AUXFUNCS_H