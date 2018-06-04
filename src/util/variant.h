#pragma once

// detect build variant and put it into VARIANT variable
#if HAVE_CUDA
#ifdef GRANULARITY_TYPE
// stringification
#define variant_xstr(s) #s
#define variant_str(s) variant_xstr(s)
#define CUDA_VARIANT "CUDA " variant_str(GRANULARITY_TYPE)
#else
#define CUDA_VARIANT "CUDA"
#endif
#endif
#if _OPENMP && HAVE_CUDA
#define VARIANT "OpenMP/" CUDA_VARIANT
#elif _OPENMP
#define VARIANT "OpenMP"
#elif HAVE_CUDA
#define VARIANT CUDA_VARIANT
#else
#define VARIANT "single-threaded CPU"
#endif
