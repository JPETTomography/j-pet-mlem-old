#pragma once

#include <cuda_runtime.h>
#include <type_traits>

#include "debug.h"

namespace util {

/// CUDA related utility classes and functions

namespace cuda {

template <typename T> inline T* device_alloc(size_t bytes) {
  T* device_ptr;
  cudaMalloc((void**)&device_ptr, bytes);
  return device_ptr;
}

/// Provides automatic management of host-device memory pointers
template <typename T> class memory {
  memory(size_t size, size_t bytes)
      : size(size),
        bytes(bytes),
        host_ptr(new T[size]),
        device_ptr(device_alloc<T>(bytes)) {}

 public:
  memory(size_t size) : memory(size, size * sizeof(T)) {}

  ~memory() {
    cudaFree(device_ptr);
    delete[] host_ptr;
  }

  void copy_to_device() {
    cudaMemcpy(device_ptr, host_ptr, bytes, cudaMemcpyHostToDevice);
  }

  void copy_from_device() {
    cudaMemcpy(host_ptr, device_ptr, bytes, cudaMemcpyDeviceToHost);
  }

  void sync_copy_from_device() {
    cudaThreadSynchronize();
    copy_from_device();
  }

  void set_on_device(int value) { cudaMemset(device_ptr, value, bytes); }

  void zero_on_device() { set_on_device(0); }

  T& operator[](size_t i) { return host_ptr[i]; }
  const T& operator[](size_t i) const { return host_ptr[i]; }

  operator T*() { return device_ptr; }

  T* const host_ptr;
  T* const device_ptr;
  const size_t size;
  const size_t bytes;
};

/// Copies automatically given data to device
template <typename T> class on_device {
  on_device(const T& on_host, size_t bytes)
      : device_ptr(device_alloc<T>(bytes)) {
    cudaMemcpy(device_ptr, &on_host, bytes, cudaMemcpyHostToDevice);
  }

 public:
  on_device(const T& on_host) : on_device(on_host, sizeof(T)) {}

  ~on_device() { cudaFree(device_ptr); }

  operator T*() { return device_ptr; }

  T* const device_ptr;
};

/// Provides underlying storage and fast copy using \c = operator

/// This can for example provide `__shared__` storage and copy data from global
/// memory with:
/// \code
/// __shared__ cuda::copy<Data> data_shared_storage;
///
/// // (1) copy Data from global memory
/// data_shared_storage = data_global_ptr;
///
/// // (2) get reference to data in shared memory
/// Data& data = *data_shared_storage;
/// \endcode
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
