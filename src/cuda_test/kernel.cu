#include <cuda.h>

static const int blocksize = 16;

__global__ void hello(char* a, int* b) { a[threadIdx.x] += b[threadIdx.x]; }

void run_kernel(char* str, int* val, int str_size, int val_size) {

  // CUDA side variables
  char* cuda_str;
  int* cuda_val;

  cudaMalloc((void**)&cuda_str, str_size);
  cudaMalloc((void**)&cuda_val, val_size);
  cudaMemcpy(cuda_str, str, str_size, cudaMemcpyHostToDevice);
  cudaMemcpy(cuda_val, val, val_size, cudaMemcpyHostToDevice);

  dim3 dimBlock(blocksize, 1);
  dim3 dimGrid(1, 1);
  hello << <dimGrid, dimBlock>>> (cuda_str, cuda_val);

  cudaMemcpy(str, cuda_str, str_size, cudaMemcpyDeviceToHost);

  cudaFree(cuda_str);
  cudaFree(cuda_str);
}
