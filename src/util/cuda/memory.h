#pragma once

#include <cuda_runtime.h>
#include <type_traits>

namespace util {
namespace cuda {

template <typename T> class copy {
  using storage_type =
      typename std::aligned_storage<sizeof(T), alignof(T)>::type;

 public:
  __device__ copy& operator=(const T* ptr) {
    const int n_blocks = (sizeof(T) + blockDim.x - 1) / blockDim.x;
    for (int block = 0; block < n_blocks; ++block) {
      const int index = blockDim.x * block + threadIdx.x;
      if (index < sizeof(T)) {
        ((char*)&storage)[index] = ((char*)ptr)[index];
      }
    }
    __syncthreads();
    return *this;
  }
  __device__ T& operator*() { return *reinterpret_cast<T*>(&storage); }

 private:
  storage_type storage;
};

}  // cuda
}  // util
