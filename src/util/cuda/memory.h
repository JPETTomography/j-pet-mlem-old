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

/// \cond PRIVATE
// Forward declarations
template <typename T> class memory;
template <typename T> class memory2D;
template <typename T> class on_device;
template <typename T> class on_device2D;
/// \endcond

/// Host-device memory.
////
/// It can be implicitly casted to device pointer. So using \c device_ptr is not
/// necessary.
template <typename T> class memory {
  memory(size_t size, size_t bytes)
      : size(size),
        bytes(bytes),
        host_ptr(new T[size]),
        device_ptr(device_alloc<T>(bytes)) {}

 protected:
  memory(size_t size, size_t bytes, T* device_ptr)
      : size(size),
        bytes(bytes),
        host_ptr(new T[size]),
        device_ptr(device_ptr) {}

 public:
  /// Allocate memory of given size on both host & device.
  memory(size_t size) : memory(size, size * sizeof(T)) {}

  ~memory() {
    cudaFree(device_ptr);
    delete[] host_ptr;
  }

  /// Copy host memory to device.
  void copy_to_device() {
    cudaMemcpy(device_ptr, host_ptr, bytes, cudaMemcpyHostToDevice);
  }

  /// Copy host memory to other on device memory.
  template <class OtherMemory> void copy_to_device(OtherMemory& other) {
    cudaMemcpy(other.device_ptr, host_ptr, bytes, cudaMemcpyHostToDevice);
  }

  /// Copy device memory to host.
  void copy_from_device() {
    cudaMemcpy(host_ptr, device_ptr, bytes, cudaMemcpyDeviceToHost);
  }

  /// Set device memory to given value.
  ////
  /// \note This sets single byte value, not float or any other type.
  void set_on_device(int value) { cudaMemset(device_ptr, value, bytes); }

  /// Zero device memory.
  void zero_on_device() { set_on_device(0); }

  T& operator[](size_t i) { return host_ptr[i]; }
  const T& operator[](size_t i) const { return host_ptr[i]; }

  T* begin() { return host_ptr; }
  T* end() { return &host_ptr[size]; }

  operator T*() { return device_ptr; }

  T* const host_ptr;
  T* const device_ptr;
  const size_t size;
  const size_t bytes;
};

/// Host-device 2D texture backed memory.
template <typename T> class memory2D : public memory<T> {
 public:
  using texture_type = texture<T, 2, cudaReadModeElementType>;

 private:
  memory2D(texture_type& tex,
           size_t width,
           size_t height,
           size_t size,
           size_t bytes,
           size_t host_pitch)
      : memory<T>(size, bytes, nullptr),
        tex(tex),
        width(width),
        height(height),
        host_pitch(host_pitch),
        device_pitch(0) {
    cudaMallocPitch((void**)&(this->device_ptr),
                    (size_t*)&device_pitch,
                    host_pitch,
                    height);
    cudaChannelFormatDesc desc = cudaCreateChannelDesc<T>();
    cudaBindTexture2D(
        NULL, &tex, this->device_ptr, &desc, width, height, device_pitch);
  }

  memory2D(texture_type& tex, size_t width, size_t height, size_t size)
      : memory2D(tex,
                 width,
                 height,
                 size,
                 size * sizeof(T),
                 width * sizeof(T)) {}

 public:
  /// Allocate memory of given dimensions on both host & device and bind to
  /// given texture.
  memory2D(texture_type& tex, size_t width, size_t height)
      : memory2D(tex, width, height, width * height) {}

  ~memory2D() { cudaUnbindTexture(&tex); }

  /// Copy from other on device memory buffer.
  memory2D& operator=(const on_device<T>& other) {
    cudaMemcpy2D(this->device_ptr,
                 device_pitch,
                 other.device_ptr,
                 host_pitch,
                 host_pitch,
                 height,
                 cudaMemcpyDeviceToDevice);
    return *this;
  }

  const texture_type& tex;
  const size_t width;
  const size_t height;
  const size_t host_pitch;
  const size_t device_pitch;
};

/// Host-device 3D texture backed memory.
template <typename T> class memory3D : public memory<T> {
 public:
  using texture_type = texture<T, 3, cudaReadModeElementType>;

 private:
  memory3D(texture_type& tex,
           size_t width,
           size_t height,
           size_t depth,
           size_t size,
           size_t bytes)
      : memory<T>(size, bytes, nullptr),
        tex(tex),
        volumeSize(make_cudaExtent(width, height, depth)) {
    cudaChannelFormatDesc desc = cudaCreateChannelDesc<T>();
    cudaMalloc3DArray((cudaArray_t*)&this->device_ptr, &desc, volumeSize);
    cudaBindTextureToArray(tex, (cudaArray_t) this->device_ptr, desc);
  }

  memory3D(texture_type& tex,
           size_t width,
           size_t height,
           size_t depth,
           size_t size)
      : memory3D(tex, width, height, depth, size, size * sizeof(T)) {}

  /// Private helper method, since 3D texture binds to special array not to
  /// linear memory, we need to handle copying differently than other classes.
  void copy_to_device(void* src_ptr, enum cudaMemcpyKind kind) {
    cudaMemcpy3DParms params = { 0 };
    params.srcPtr = make_cudaPitchedPtr(src_ptr,
                                        volumeSize.width * sizeof(T),
                                        volumeSize.width,
                                        volumeSize.height);
    params.dstArray = (cudaArray_t) this->device_ptr;
    params.extent = volumeSize;
    params.kind = kind;
    cudaMemcpy3D(&params);
  }

 public:
  /// Allocate memory of given dimensions on both host & device and bind to
  /// given texture.
  memory3D(texture_type& tex, size_t width, size_t height, size_t depth)
      : memory3D(tex, width, height, depth, width * height * depth) {}

  ~memory3D() { cudaUnbindTexture(&tex); }

  /// Copy host memory to device.
  void copy_to_device() {
    copy_to_device(this->host_ptr, cudaMemcpyHostToDevice);
  }

  /// Copy from other on device memory buffer.
  memory3D& operator=(const on_device<T>& other) {
    copy_to_device(other.device_ptr, cudaMemcpyDeviceToDevice);
    return *this;
  }

  /// Copy device memory to host is NOT SUPPORTED for 3D memory.
  void copy_from_device() {}

  const texture_type& tex;
  const cudaExtent volumeSize;
};

/// Device only memory.
////
/// It can be implicitly casted to device pointer. So using \c device_ptr is not
/// necessary.
template <typename T> class on_device {
  on_device(const T& on_host, size_t bytes)
      : bytes(bytes), device_ptr(device_alloc<T>(bytes)) {
    cudaMemcpy(device_ptr, &on_host, bytes, cudaMemcpyHostToDevice);
  }

  on_device(const T* on_host, size_t, size_t bytes)
      : bytes(bytes), device_ptr(device_alloc<T>(bytes)) {
    cudaMemcpy(device_ptr, on_host, bytes, cudaMemcpyHostToDevice);
  }

  on_device(size_t, size_t bytes)
      : bytes(bytes), device_ptr(device_alloc<T>(bytes)) {}

 protected:
  on_device(size_t bytes, T* device_ptr)
      : bytes(bytes), device_ptr(device_ptr) {}

 public:
  /// Allocate and copy from host.
  on_device(const T& on_host) : on_device(on_host, sizeof(T)) {}

  /// Allocate and copy as array of given size from host.
  on_device(const T* on_host, size_t size)
      : on_device(on_host, size, sizeof(T) * size) {}

  /// Just allocate on device.
  on_device() : bytes(sizeof(T)), device_ptr(device_alloc<T>(sizeof(T))) {}

  /// Just allocate as array of given size on device.
  on_device(size_t size) : on_device(size, size * sizeof(T)) {}

  ~on_device() { cudaFree(device_ptr); }

  /// Set device memory to given value.
  void set_on_device(int value) { cudaMemset(device_ptr, value, bytes); }

  /// Zero device memory.
  void zero_on_device() { set_on_device(0); }

  operator T*() { return device_ptr; }

  T* const device_ptr;
  const size_t bytes;
};

/// Device only 2D texture backed memory.
template <typename T> class on_device2D : public on_device<T> {
 public:
  using texture_type = texture<T, 2, cudaReadModeElementType>;

 private:
  on_device2D(texture_type& tex,
              size_t width,
              size_t height,
              size_t bytes,
              size_t host_pitch)
      : on_device<T>(bytes, nullptr), tex(tex) {
    size_t device_pitch;
    cudaMallocPitch((void**)&(this->device_ptr),
                    (size_t*)&device_pitch,
                    host_pitch,
                    height);
    cudaChannelFormatDesc desc = cudaCreateChannelDesc<T>();
    cudaBindTexture2D(
        NULL, &tex, this->device_ptr, &desc, width, height, device_pitch);
  }

  on_device2D(texture_type& tex, size_t width, size_t height, size_t size)
      : on_device2D(tex, width, height, size * sizeof(T), width * sizeof(T)) {}

  on_device2D(texture_type& tex,
              T* on_host,
              size_t width,
              size_t height,
              size_t bytes,
              size_t host_pitch)
      : on_device<T>(bytes, nullptr), tex(tex) {
    size_t device_pitch;
    cudaMallocPitch((void**)&(this->device_ptr),
                    (size_t*)&device_pitch,
                    host_pitch,
                    height);
    cudaMemcpy2D(this->device_ptr,
                 device_pitch,
                 on_host,
                 host_pitch,
                 host_pitch,
                 height,
                 cudaMemcpyHostToDevice);
    cudaChannelFormatDesc desc = cudaCreateChannelDesc<T>();
    cudaBindTexture2D(
        NULL, &tex, this->device_ptr, &desc, width, height, device_pitch);
  }

  on_device2D(texture_type& tex,
              T* on_host,
              size_t width,
              size_t height,
              size_t size)
      : on_device2D(tex,
                    on_host,
                    width,
                    height,
                    size * sizeof(T),
                    width * sizeof(T)) {}

 public:
  /// Just allocate with given dimensions on device and bind to texture.
  on_device2D(texture_type& tex, size_t width, size_t height)
      : on_device2D(tex, width, height, width * height) {}

  /// Allocate with given dimensions, copy from host and bind to texture.
  on_device2D(texture_type& tex, T* on_host, size_t width, size_t height)
      : on_device2D(tex, on_host, width, height, width * height) {}

  ~on_device2D() { cudaUnbindTexture(&tex); }

  const texture_type& tex;
};

/// Provides underlying storage and fast copy using \c = operator
////
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
