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
void populateDiagonalThomasAlgorithm(real *la, real *lb, real *lc, int N, real phi)
{
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

// TODO: Correct this function
void prefactorizeThomas(real *la, real *lb, real *lc, real *c_prime, real *denominator, int N)
{
    c_prime[0] = lc[0] / lb[0];
    denominator[0] = 1.0f / lb[0];

    for (int i = 1; i < N; i++)
    {
        real denom = 1.0f / (lb[i] - c_prime[i - 1] * la[i]);
        if (i < N - 1)
        {
            c_prime[i] = lc[i] * denom;
        }
        denominator[i] = denom;
    }
}

void createDirectories(real delta_t, real delta_x, real delta_y, char *pathToSaveData)
{
    // Build the path
    char path[MAX_STRING_SIZE];

#ifndef CABLEEQ

    snprintf(path, MAX_STRING_SIZE, "./simulation_files/dt_%.5g_dx_%.5g_dy_%.5g/%s/%s/%s/%s/%s", delta_t, delta_x, delta_y, EXECUTION_TYPE, REAL_TYPE, PROBLEM, CELL_MODEL, METHOD);

#else // if def CABLEEQ

    snprintf(path, MAX_STRING_SIZE, "./simulation_files/dt_%.5g_dx_%.5g/%s/%s/%s/%s/%s", delta_t, delta_x, EXECUTION_TYPE, REAL_TYPE, PROBLEM, CELL_MODEL, METHOD);

#endif // not CABLEEQ

#ifdef THETA

    // Add theta to the path
    char thetaPath[MAX_STRING_SIZE];
    snprintf(thetaPath, MAX_STRING_SIZE, "%.2lf", THETA);
    strcat(path, "/");
    strcat(path, thetaPath);

#endif // THETA

    // Update pathToSaveData
    strcpy(pathToSaveData, path);

    // Create directories
    char temp[MAX_STRING_SIZE];
    strcpy(temp, path);

    char *p = temp;
    while (*p) {
        if (*p == '/') {
            *p = '\0';
            if (mkdir(temp, 0777) != 0 && errno != EEXIST)
                ERRORMSG("Error creating dir %s: %s\n", temp, strerror(errno));

            *p = '/';
        }
        p++;
    }
    if (mkdir(temp, 0777) != 0 && errno != EEXIST)
        ERRORMSG("Error creating dir %s: %s\n", temp, strerror(errno));
}

void initializeTimeArray(real *timeArray, int M, real dt)
{
    for (int i = 0; i < M; i++)
        timeArray[i] = i * dt;
}

int lim(int num, int N)
{
    if (num == -1)
        return 1;
    else if (num == N)
        return N - 2;
    return num;
}

#ifdef MV

real rescaleVm(real Vm)
{
    return 85.7f * Vm - 84.0f;
}

#endif // MV

#ifndef CONVERGENCE_ANALYSIS_FORCING_TERM
#if defined(MONODOMAIN) || defined(CABLEEQ)

void populateStimuli(Stimulus *stimuli, real delta_x, real delta_y)
{
    for (int i = 0; i < numberOfStimuli; i++)
    {
        stimuli[i].amplitude = stimuliAmplitude;
        stimuli[i].begin = stimuliBegin[i];
        stimuli[i].duration = stimuliDuration;

        // Discretized limits of stimulation areas
        stimuli[i].xMaxDisc = round(stimulixMax[i] / delta_x);
        stimuli[i].xMinDisc = round(stimulixMin[i] / delta_x);
        stimuli[i].yMaxDisc = round(stimuliyMax[i] / delta_y);
        stimuli[i].yMinDisc = round(stimuliyMin[i] / delta_y);

#ifdef CABLEEQ

        // Only one stimulus for CABLEEQ
        if (i > 0)
            stimuli[i].amplitude = 0.0f;

#endif // CABLEEQ
    }
}

#endif // MONODOMAIN || CABLEEQ
#endif // not CONVERGENCE_ANALYSIS_FORCING_TERM

#if defined(SERIAL) || defined(OPENMP)

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

real calculateNorm2Error(real **Vm, real **exact, int Nx, int Ny, real totalTime, real delta_x, real delta_y)
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
            sum += ((Vm[i][j] - solution) * (Vm[i][j] - solution));
        }
    }
    return sqrt(sum / (Nx * Ny));
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

void saveFrame(FILE *file, real **Vm, int Nx, int Ny)
{
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
#ifndef MV
            fprintf(file, "%e ", Vm[i][j]);
#else  // if def MV
            fprintf(file, "%e ", rescaleVm(Vm[i][j]));
#endif // not MV
        }
        fprintf(file, "\n");
    }
}

#else // if def CABLEEQ

void saveFrame(FILE *file, real *Vm, int N)
{
    for (int i = 0; i < N; i++)
    {
#ifndef MV
        fprintf(file, "%e ", Vm[i]);
#else  // if def MV
        fprintf(file, "%e ", rescaleVm(Vm[i]));
#endif // not MV
    }
    fprintf(file, "\n");
}

#endif // not CABLEEQ
#endif // SERIAL || OPENMP

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

void shift2DVariableToLeft(real *Var, int Nx, int Ny, real length, real delta_x, real delta_y, real initValue, char *varName)
{
    real *temp = (real *)malloc(Nx * sizeof(real));
    int index;
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            index = i * Nx + j;
            temp[j] = Var[index];
        }

        int lengthIndex = round(length / delta_x) + 1;
        for (int j = 0; j < Nx - lengthIndex; j++)
        {
            index = i * Nx + j;
            Var[index] = temp[j + lengthIndex];
        }
        for (int j = Nx - lengthIndex; j < Nx; j++)
        {
            index = i * Nx + j;
            Var[index] = initValue;
        }
    }
    free(temp);
    SUCCESSMSG("Variable %s shifted to the left by %.2f cm\n", varName, length);
}

void saveFrame(FILE *file, real *Vm, int Nx, int Ny)
{
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            int index = i * Nx + j;
#ifndef MV
            fprintf(file, "%e ", Vm[index]);
#else  // if def MV
            fprintf(file, "%e ", rescaleVm(Vm[index]));
#endif // not MV
        }
        fprintf(file, "\n");
    }
}

#endif // GPU

#endif // AUXFUNCS_H