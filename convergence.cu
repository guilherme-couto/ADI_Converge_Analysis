#include "./methods.h"

int main(int argc, char *argv[])
{
    // Parameters
    char *method;
    real delta_t;
    real delta_x;
    real theta;

    // Read parameters
    if (argc != 5)
    {
        printf("Usage: %s method delta_t delta_x theta\n", argv[0]);
        return 1;
    }
    method = argv[1];
    delta_t = atof(argv[2]);
    delta_x = atof(argv[3]);
    theta = atof(argv[4]);

    int N = 3;

    real **RHS = (real **)malloc(N * sizeof(real *));
    real *c_prime = (real *)malloc(N * sizeof(real));   // aux for Thomas
    real *d_prime = (real *)malloc(N * sizeof(real));   // aux for Thomas
    for (int i = 0; i < N; i++)
    {
        RHS[i] = (real *)malloc(N * sizeof(real));
    }

    // Populate RHS
    RHS[0][0] = 13;
    RHS[0][1] = 14;
    RHS[0][2] = 15;
    RHS[1][0] = 4;
    RHS[1][1] = 5;
    RHS[1][2] = 6;
    RHS[2][0] = -5;
    RHS[2][1] = -4;
    RHS[2][2] = -3;

    printf("\nMatriz RHS\n");
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            printf("%lf ", RHS[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    // Auxiliary arrays for Thomas algorithm
    real *la = (real *)malloc(N * sizeof(real));
    real *lb = (real *)malloc(N * sizeof(real));
    real *lc = (real *)malloc(N * sizeof(real));
    real *d = (real *)malloc(N * sizeof(real));

    populateDiagonalThomasAlgorithm(la, lb, lc, N, 2);

    for (int i = 0; i < N; i++)
        thomasAlgorithm(la, lb, lc, c_prime, d_prime, N, RHS[i]);

    printf("\nResultado do Thomas para diff imp x\n");
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            printf("%lf ", RHS[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    for (int j = 0; j < N; j++)
    {
        copyColumnToVector(RHS, d, N, j);
        thomasAlgorithm(la, lb, lc, c_prime, d_prime, N, d);
        copyVectorToColumn(RHS, d, N, j);
    }

    printf("\nResultado do Thomas para diff imp y\n");
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            printf("%lf ", RHS[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    // Free memory from host
    for (int i = 0; i < N; i++)
        free(RHS[i]);
    free(RHS);
    free(la);
    free(lb);
    free(lc);
    free(c_prime);
    free(d_prime);
    free(d);

    return 0;

    
    // Call function
    runSimulation(method, delta_t, delta_x, theta);

    return 0;
}