#pragma once

// detect build variant and put it into VARIANT variable
#if _OPENMP && HAVE_CUDA
#define VARIANT "OpenMP/CUDA"
#elif _OPENMP
#define VARIANT "OpenMP"
#elif HAVE_CUDA
#define VARIANT "CUDA"
#else
#define VARIANT "single-threaded CPU"
#endif
