// CUDA functions debug automatic wrapper
//
// Author:
//   Adam Strzelecki <adam.strzelecki@uj.edu.pl>
//
// Discussion:
//   This header wraps commonly used functions (listed below) that may generate
//   errors into conditional checking for result printing error message and
//   aborting execution when error code is returned.

#pragma once

static cudaError cudbgLastError;

#define CUDBG(f, ...)                                           \
  if ((cudbgLastError = cuda##f(__VA_ARGS__)) != cudaSuccess) { \
    fprintf(stderr,                                             \
            "%s:%d cuda%s() %d: %s\n",                          \
            __FILE__,                                           \
            __LINE__,                                           \
            #f,                                                 \
            cudbgLastError,                                     \
            cudaGetErrorString(cudbgLastError));                \
    abort();                                                    \
  }

// automatically wraps all commonly used CUDA functions with error check

#if __CUDACC__

#define cudaBindSurfaceToArray(...) CUDBG(BindSurfaceToArray, __VA_ARGS__)
#define cudaCreateTextureObject(...) CUDBG(CreateTextureObject, __VA_ARGS__)
#define cudaDestroyTextureObject(...) CUDBG(DestroyTextureObject, __VA_ARGS__)
#define cudaFree(...) /* -------------> */ CUDBG(Free, __VA_ARGS__)
#define cudaMalloc(...) /* -----------> */ CUDBG(Malloc, __VA_ARGS__)
#define cudaMallocArray(...) /* ------> */ CUDBG(MallocArray, __VA_ARGS__)
#define cudaMallocPitch(...) /* ------> */ CUDBG(MallocPitch, __VA_ARGS__)
#define cudaMemcpy(...) /* -----------> */ CUDBG(Memcpy, __VA_ARGS__)
#define cudaMemcpyToArray(...) /* ----> */ CUDBG(MemcpyToArray, __VA_ARGS__)
#define cudaMemcpy2D(...) /* ---------> */ CUDBG(Memcpy2D, __VA_ARGS__)
#define cudaMemset(...) /* -----------> */ CUDBG(Memset, __VA_ARGS__)
#define cudaSetDevice(...) /* --------> */ CUDBG(SetDevice, __VA_ARGS__)
#define cudaThreadSynchronize(...) /* > */ CUDBG(ThreadSynchronize, __VA_ARGS__)

#endif
