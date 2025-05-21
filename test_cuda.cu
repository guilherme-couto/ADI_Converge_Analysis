#include <cuda_runtime.h>
#include <stdio.h>

int main() {
    cudaDeviceProp prop;
    int deviceCount;
    
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount > 0) {
        cudaGetDeviceProperties(&prop, 0);
        printf("Device name: %s\n", prop.name);
    } else {
        printf("No CUDA devices found.\n");
    }

    return 0;
}
