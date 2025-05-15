#ifndef AUXFUNCS_H
#define AUXFUNCS_H

#include "core_definitions.h"
#include "config_parser.h"
#include "logger.h"

void populateDiagonalThomasAlgorithm(real *la, real *lb, real *lc, int N, real phi);
void prefactorizeThomas(real *la, real *lb, real *lc, real *c_prime, real *denominator, int N); // TODO: Correct this function
const int createDirectories(char *dir_path, bool remove_old_files);
void initializeTimeArray(real *timeArray, int M, real dt);
void initializeMeasurement(Measurement *measurement);
const int lim(int num, int N);
const real rescaleVm(real Vm);
int populateStimuli(SimulationConfig *config);
void saveCopyOfSimulationConfig(const char *ini_file_path, const char *output_dir);
void tridiagonalSystemSolver(real *la, real *lb, real *lc, real *c_prime, real *d_prime, int N, real *d, real *result);
const int saveSimulationInfos(const SimulationConfig *config, const Measurement *measurement);

real calculateNorm2Error(real *Vm, real *exact, int Nx, int Ny, real totalTime, real delta_x, real delta_y, real Lx, real Ly);
void initializeVariableWithExactSolution(real *Var, int Nx, int Ny, real delta_x, real delta_y, real Lx, real Ly);
void initializeVariableWithValue(real *Var, int Nx, int Ny, real value);
void initializeVariableFromFile(real *Var, int Nx, int Ny, char *filename, real delta_x, real delta_y, char *varName, real reference_dx, real reference_dy, real Lx, real Ly);
void shiftVariableToLeft(real *Var, int Nx, int Ny, real length, real delta_x, real delta_y, real initValue, char *varName);

#ifdef CABLEEQ

void initialize1DVariableWithValue(real *Var, int N, real value);
void initialize1DVariableFromFile(real *Var, int N, char *filename, real delta_x, char *varName, real reference_dx, real Lx);
void shift1DVariableToLeft(real *Var, int N, real length, real delta_x, real initValue, char *varName);

#endif // CABLEEQ

#endif // AUXFUNCS_H