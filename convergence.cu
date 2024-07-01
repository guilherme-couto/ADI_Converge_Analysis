#include "./src/methods.h"

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

    // Call function
    runSimulation(method, delta_t, delta_x, theta);

    return 0;

    // TESTES THOMAS ///////////////////////////////////////////////////////////////////////////////////////////////////////
    // int N = 3;

    // real **RHS = (real **)malloc(N * sizeof(real *));
    // real *c_prime = (real *)malloc(N * sizeof(real));   // aux for Thomas
    // real *d_prime = (real *)malloc(N * sizeof(real));   // aux for Thomas
    // for (int i = 0; i < N; i++)
    // {
    //     RHS[i] = (real *)malloc(N * sizeof(real));
    // }

    // // Populate RHS
    // RHS[0][0] = 13;
    // RHS[0][1] = 14;
    // RHS[0][2] = 15;
    // RHS[1][0] = 4;
    // RHS[1][1] = 5;
    // RHS[1][2] = 6;
    // RHS[2][0] = -5;
    // RHS[2][1] = -4;
    // RHS[2][2] = -3;

    // printf("\nMatriz RHS\n");
    // for (int i = 0; i < N; i++)
    // {
    //     for (int j = 0; j < N; j++)
    //     {
    //         printf("%lf ", RHS[i][j]);
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    // // Auxiliary arrays for Thomas algorithm
    // real *la = (real *)malloc(N * sizeof(real));
    // real *lb = (real *)malloc(N * sizeof(real));
    // real *lc = (real *)malloc(N * sizeof(real));
    // real *d = (real *)malloc(N * sizeof(real));
    // real *solution = (real *)malloc(N * sizeof(real));

    // populateDiagonalThomasAlgorithm(la, lb, lc, N, 2);

    // for (int i = 0; i < N; i++)
    //     thomasAlgorithm(la, lb, lc, c_prime, d_prime, N, RHS[i]);

    // printf("\nResultado do Thomas para diff imp x\n");
    // for (int i = 0; i < N; i++)
    // {
    //     for (int j = 0; j < N; j++)
    //     {
    //         printf("%lf ", RHS[i][j]);
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    // for (int j = 0; j < N; j++)
    // {
    //     copyColumnToVector(RHS, d, N, j);
    //     thomasAlgorithm(la, lb, lc, c_prime, d_prime, N, d);
    //     copyVectorToColumn(RHS, d, N, j);
    // }

    // printf("\nResultado do Thomas para diff imp y\n");
    // for (int i = 0; i < N; i++)
    // {
    //     for (int j = 0; j < N; j++)
    //     {
    //         printf("%lf ", RHS[i][j]);
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    // // Free memory from host
    // for (int i = 0; i < N; i++)
    //     free(RHS[i]);
    // free(RHS);
    // free(la);
    // free(lb);
    // free(lc);
    // free(c_prime);
    // free(d_prime);
    // free(d);
    // free(solution);

    // return 0;
    // /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // TESTE PARALELO /////////////////////////////////////////////////////////////////////////////////////////
    // int N = 3;

    // real *RHS;
    // RHS = (real *)malloc(N * N * sizeof(real));
    // RHS[0]=13;
    // RHS[1]=14;
    // RHS[2]=15;
    // RHS[3]=4;
    // RHS[4]=5;
    // RHS[5]=6;
    // RHS[6]=-5;
    // RHS[7]=-4;
    // RHS[8]=-3;

    // // Auxiliary arrays for Thomas algorithm
    // real *la = (real *)malloc(N * sizeof(real));
    // real *lb = (real *)malloc(N * sizeof(real));
    // real *lc = (real *)malloc(N * sizeof(real));

    
    // populateDiagonalThomasAlgorithm(la, lb, lc, N, 2);

    // // Prefactorization
    // thomasFactorConstantBatch(la, lb, lc, N);

    // real *d_RHS, *d_aux;
    // real *d_la, *d_lb, *d_lc;
    // cudaError_t cudaStatus1, cudaStatus2, cudaStatus3, cudaStatus6, cudaStatus7, cudaStatus8;
    
    // cudaStatus2 = cudaMalloc(&d_RHS, N * N * sizeof(real));
    // cudaStatus2 = cudaMalloc(&d_aux, N * N * sizeof(real));
    // cudaStatus6 = cudaMalloc(&d_la, N * sizeof(real));
    // cudaStatus7 = cudaMalloc(&d_lb, N * sizeof(real));
    // cudaStatus8 = cudaMalloc(&d_lc, N * sizeof(real));
    // if (cudaStatus2 != cudaSuccess || cudaStatus6 != cudaSuccess || cudaStatus7 != cudaSuccess || cudaStatus8 != cudaSuccess)
    // {
    //     printf("cudaMalloc failed!\n");
    //     exit(EXIT_FAILURE);
    // }

    // // Copy memory from host to device of the matrices (2D arrays)
    // cudaStatus1 = cudaMemcpy(d_RHS, RHS, N * N * sizeof(real), cudaMemcpyHostToDevice);
    // if (cudaStatus1 != cudaSuccess)
    // {
    //     printf("cudaMemcpy failed!\n");
    //     exit(EXIT_FAILURE);
    // }

    // // Copy memory of diagonals from host to device
    // cudaStatus1 = cudaMemcpy(d_la, la, N * sizeof(real), cudaMemcpyHostToDevice);
    // cudaStatus2 = cudaMemcpy(d_lb, lb, N * sizeof(real), cudaMemcpyHostToDevice);
    // cudaStatus3 = cudaMemcpy(d_lc, lc, N * sizeof(real), cudaMemcpyHostToDevice);
    // if (cudaStatus1 != cudaSuccess || cudaStatus2 != cudaSuccess || cudaStatus3 != cudaSuccess)
    // {
    //     printf("cudaMemcpy failed!\n");
    //     exit(EXIT_FAILURE);
    // }
    // printf("All cudaMallocs done!\n");

    // // Block and grid size
    // // For parallel Thomas
    // printf("N = %d\n", N);
    // int numBlocks = N / 100;
    // if (numBlocks == 0)
    //     numBlocks = 1;
    // int blockSize = round(N / numBlocks) + 1;
        
    // if (blockSize % 32 != 0)
    //     blockSize = 32 * ((blockSize / 32) + 1);

    // // For other kernels
    // int GRID_SIZE = ceil((N*N*1.0) / (BLOCK_SIZE*1.0));
    // if (GRID_SIZE == 0)
    //     GRID_SIZE = 1;

    // printf("For 1st Part and Transpose -> Grid size %d, Block size %d\n", GRID_SIZE, BLOCK_SIZE);
    // printf("Total for 1st Part and Transpose: %d\n", GRID_SIZE*BLOCK_SIZE);
    // printf("For Thomas Algorithm -> Grid size: %d, Block size: %d\n", numBlocks, blockSize);
    // printf("Total for Thomas Algorithm: %d\n", numBlocks*blockSize);

    // cuThomasConstantBatch<<<numBlocks, blockSize>>>(d_la, d_lb, d_lc, d_RHS, N);
    // cudaDeviceSynchronize();

    // //Copy d_RHS to RHS
    // cudaStatus1 = cudaMemcpy(RHS, d_RHS, N * N * sizeof(real), cudaMemcpyDeviceToHost);
    // if (cudaStatus1 != cudaSuccess)
    // {
    //     printf("cudaMemcpy failed!\n");
    //     exit(EXIT_FAILURE);
    // }

    // printf("\nResultado sem transpor a matriz RHS\n");
    // for (int i = 0; i < N; i++)
    // {
    //     for (int j = 0; j < N; j++)
    //     {
    //         int index = i * N + j;
    //         printf("%lf ", RHS[index]);
    //     }
    //     printf("\n");
    // }

    // //########################
    // RHS[0]=13;
    // RHS[1]=14;
    // RHS[2]=15;
    // RHS[3]=4;
    // RHS[4]=5;
    // RHS[5]=6;
    // RHS[6]=-5;
    // RHS[7]=-4;
    // RHS[8]=-3;

    // cudaStatus1 = cudaMemcpy(d_RHS, RHS, N * N * sizeof(real), cudaMemcpyHostToDevice);

    // transposeDiagonalCol<<<GRID_SIZE, BLOCK_SIZE>>>(d_RHS, d_aux, N);
    // cudaDeviceSynchronize();

    // //Copy d_RHS to RHS
    // cudaStatus1 = cudaMemcpy(RHS, d_aux, N * N * sizeof(real), cudaMemcpyDeviceToHost);
    
    // printf("\nRHS transposta (interleaved)\n");
    // for (int i = 0; i < N; i++)
    // {
    //     for (int j = 0; j < N; j++)
    //     {
    //         int index = i * N + j;
    //         printf("%lf ", RHS[index]);
    //     }
    //     printf("\n");
    // }

    // cuThomasConstantBatch<<<numBlocks, blockSize>>>(d_la, d_lb, d_lc, d_aux, N);
    // cudaDeviceSynchronize();

    // //Copy d_RHS to RHS
    // cudaStatus1 = cudaMemcpy(RHS, d_aux, N * N * sizeof(real), cudaMemcpyDeviceToHost);
    // if (cudaStatus1 != cudaSuccess)
    // {
    //     printf("cudaMemcpy failed!\n");
    //     exit(EXIT_FAILURE);
    // }

    // printf("\nResultado transpondo a matriz RHS\n");
    // for (int i = 0; i < N; i++)
    // {
    //     for (int j = 0; j < N; j++)
    //     {
    //         int index = i * N + j;
    //         printf("%lf ", RHS[index]);
    //     }
    //     printf("\n");
    // }

    // transposeDiagonalCol<<<GRID_SIZE, BLOCK_SIZE>>>(d_aux, d_RHS, N);
    // cudaDeviceSynchronize();

    // //Copy d_RHS to RHS
    // cudaStatus1 = cudaMemcpy(RHS, d_RHS, N * N * sizeof(real), cudaMemcpyDeviceToHost);
    
    // printf("\n2a RHS transposta (interleaved)\n");
    // for (int i = 0; i < N; i++)
    // {
    //     for (int j = 0; j < N; j++)
    //     {
    //         int index = i * N + j;
    //         printf("%lf ", RHS[index]);
    //     }
    //     printf("\n");
    // }

    // cuThomasConstantBatch<<<numBlocks, blockSize>>>(d_la, d_lb, d_lc, d_RHS, N);
    // cudaDeviceSynchronize();

    // //Copy d_RHS to RHS
    // cudaStatus1 = cudaMemcpy(RHS, d_RHS, N * N * sizeof(real), cudaMemcpyDeviceToHost);
    // if (cudaStatus1 != cudaSuccess)
    // {
    //     printf("cudaMemcpy failed!\n");
    //     exit(EXIT_FAILURE);
    // }

    // printf("\nResultado Pós segundo Sistema\n");
    // for (int i = 0; i < N; i++)
    // {
    //     for (int j = 0; j < N; j++)
    //     {
    //         int index = i * N + j;
    //         printf("%lf ", RHS[index]);
    //     }
    //     printf("\n");
    // }



    // // Free memory from host
    // free(RHS);
    // free(la);
    // free(lb);
    // free(lc);

    // // Free memory from device
    // cudaFree(d_RHS);
    // cudaFree(d_aux);
    // cudaFree(d_la);
    // cudaFree(d_lb);
    // cudaFree(d_lc);
    // return 0;
    
}