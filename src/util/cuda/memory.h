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

template <typename T>
inline cudaArray_t device_alloc_array(
    const size_t width,
    const size_t height,
    const cudaChannelFormatDesc desc = cudaCreateChannelDesc<T>()) {
  cudaArray_t array;
  cudaMallocArray(&array, &desc, width, height);
  return array;
}

template <typename T>
inline cudaArray_t device_alloc_array3D(
    const cudaExtent volumeSize,
    const cudaChannelFormatDesc desc = cudaCreateChannelDesc<T>()) {
  cudaArray_t array;
  cudaMalloc3DArray(&array, &desc, volumeSize);
  return array;
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
    if (device_ptr)
      cudaFree(device_ptr);
    if (host_ptr)
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

  /// Copy from other on device memory buffer.
  on_device& operator=(const memory<T>& other) {
    cudaMemcpy(device_ptr, other.device_ptr, bytes, cudaMemcpyDeviceToDevice);
    return *this;
  }

  /// Copy from other on device memory buffer.
  on_device& operator=(const on_device& other) {
    cudaMemcpy(device_ptr, other.device_ptr, bytes, cudaMemcpyDeviceToDevice);
    return *this;
  }

  operator T*() { return device_ptr; }

  T* const device_ptr;
  const size_t bytes;
};

/// Device 2D texture backed memory.
template <typename T> class texture2D {
 public:
  using texture_type = texture<T, 2, cudaReadModeElementType>;

 private:
  texture2D(texture_type& tex,
            size_t width,
            size_t height,
            cudaChannelFormatDesc desc,
            const T* values)
      : tex(tex),
        width(width),
        height(height),
        array(device_alloc_array<T>(width, height, desc)) {
    if (values) {
      copy_to_array(values, cudaMemcpyHostToDevice);
    }
    cudaBindTextureToArray(tex, array, desc);
  }

  /// Private helper method, since 3D texture binds to special array not to
  /// linear memory, we need to handle copying differently than other classes.
  void copy_to_array(const void* src_ptr, enum cudaMemcpyKind kind) {
    cudaMemcpy2DToArray(array,
                        0,
                        0,
                        src_ptr,
                        width * sizeof(T),
                        width * sizeof(T),
                        height,
                        kind);
  }

 public:
  /// Allocate memory of given dimensions on device and bind to given texture.
  texture2D(texture_type& tex,         ///< texture reference to bind to
            size_t width,              ///< width of the texture
            size_t height,             ///< height of the texture
            const T* values = nullptr  ///< host values to copy to texture
                                       ///  (pass \c nullptr to skip)
            )
      : texture2D(tex, width, height, cudaCreateChannelDesc<T>(), values) {}

  ~texture2D() {
    cudaUnbindTexture(&tex);
    cudaFreeArray(array);
  }

  /// Copy from host memory buffer.
  texture2D& operator=(const T* other_ptr) {
    copy_to_array(other_ptr, cudaMemcpyHostToDevice);
    return *this;
  }

  /// Copy from other on device memory buffer.
  texture2D& operator=(const memory<T>& other) {
    copy_to_array(other.device_ptr, cudaMemcpyDeviceToDevice);
    return *this;
  }

  /// Copy from other on device memory buffer.
  texture2D& operator=(const on_device<T>& other) {
    copy_to_array(other.device_ptr, cudaMemcpyDeviceToDevice);
    return *this;
  }

  const texture_type& tex;
  const size_t width;
  const size_t height;
  const cudaArray_t array;
};

/// Device 3D texture backed memory.
template <typename T> class texture3D {
 public:
  using texture_type = texture<T, 3, cudaReadModeElementType>;

 private:
  texture3D(texture_type& tex,
            cudaExtent volumeSize,
            cudaChannelFormatDesc desc,
            const T* values)
      : tex(tex),
        volumeSize(volumeSize),
        array(device_alloc_array3D<T>(volumeSize, desc)) {
    if (values) {
      copy_to_array(values, cudaMemcpyHostToDevice);
    }
    cudaBindTextureToArray(tex, array, desc);
  }

  /// Private helper method, since 3D texture binds to special array not to
  /// linear memory, we need to handle copying differently than other classes.
  void copy_to_array(const void* src_ptr, enum cudaMemcpyKind kind) {
    cudaMemcpy3DParms params = { 0 };
    // FIXME: for some reason make_cudaPitchedPtr wants void* not const void*
    params.srcPtr = make_cudaPitchedPtr(const_cast<void*>(src_ptr),
                                        volumeSize.width * sizeof(T),
                                        volumeSize.width,
                                        volumeSize.height);
    params.dstArray = array;
    params.extent = volumeSize;
    params.kind = kind;
    cudaMemcpy3D(&params);
  }

 public:
  /// Allocate memory of given dimensions on device and bind to given texture.
  texture3D(texture_type& tex,         ///< texture reference to bind to
            size_t width,              ///< width of the texture
            size_t height,             ///< height of the texture
            size_t depth,              ///< depth of the texture
            const T* values = nullptr  ///< host values to copy to texture
                                       ///  (pass \c nullptr to skip)
            )
      : texture3D(tex,
                  make_cudaExtent(width, height, depth),
                  cudaCreateChannelDesc<T>(),
                  values) {}

  ~texture3D() {
    cudaUnbindTexture(&tex);
    cudaFreeArray(array);
  }

  /// Copy from host memory buffer.
  texture3D& operator=(const T* other_ptr) {
    copy_to_array(other_ptr, cudaMemcpyHostToDevice);
    return *this;
  }

  /// Copy from other on device memory buffer.
  texture3D& operator=(const memory<T>& other) {
    copy_to_array(other.device_ptr, cudaMemcpyDeviceToDevice);
    return *this;
  }

  /// Copy from other on device memory buffer.
  texture3D& operator=(const on_device<T>& other) {
    copy_to_array(other.device_ptr, cudaMemcpyDeviceToDevice);
    return *this;
  }

  const texture_type& tex;
  const cudaExtent volumeSize;
  const cudaArray_t array;
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
