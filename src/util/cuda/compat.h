// CUDA generic code compatibility macros & functions
//
// Author:
//   Adam Strzelecki <adam.strzelecki@uj.edu.pl>
//
// Discussion:
//   Purpose of this header is to provide maximum share-ability between generic
//   CPU implementation and CUDA implementation. This is done using following:
//
//   (1) All mathematical function are exposed via compat:: namespace i.e.
//       compat::cos which map to stdlib or CUDA depending on build context.
//
//   (2) Special underscore macro _ is used to mark functions and methods
//       compatible with CUDA and it is replaced with __device__ __host__ when
//       building using CUDA.

#pragma once

#if !__CUDACC__
#include <cmath>
#else
#include <limits.h>
#include <float.h>
#endif

namespace compat {

#if __CUDACC__

#define _ __device__ __host__
#define constexpr

_ int min(const int a, const int b) { return ::min(a, b); }
_ int max(const int a, const int b) { return ::max(a, b); }
_ float min(const float a, const float b) { return fminf(a, b); }
_ float max(const float a, const float b) { return fmaxf(a, b); }
_ float round(const float a) { return roundf(a); }
_ float ceil(const float a) { return ceilf(a); }
_ float floor(const float a) { return floorf(a); }
_ float sqrt(const float a) { return sqrtf(a); }
_ float sin(const float a) { return sinf(a); }
_ float cos(const float a) { return cosf(a); }
_ float tan(const float a) { return tanf(a); }
_ float atan(const float a) { return atanf(a); }
_ float atan2(const float a, const float b) { return atan2f(a, b); }
_ float pow(const float a, const float b) { return powf(a, b); }
_ float exp(const float a) { return expf(a); }
_ float log(const float a) { return logf(a); }

template <typename F> _ constexpr F numeric_max();
template <> _ constexpr float numeric_max<float>() { return FLT_MAX; }
template <> _ constexpr int numeric_max<int>() { return INT_MAX; }
template <> _ constexpr unsigned int numeric_max<unsigned int>() {
  return UINT_MAX;
}
template <> _ constexpr unsigned long numeric_max<unsigned long>() {
  return ULONG_MAX;
}

#else

#define _

#if _MSC_VER
template <typename F> F min(const F a, const F b) { return (a < b) ? a : b; }
template <typename F> F max(const F a, const F b) { return (a > b) ? a : b; }
#else
template <typename F> F min(const F a, const F b) { return std::min(a, b); }
template <typename F> F max(const F a, const F b) { return std::max(a, b); }
#endif
template <typename F> F round(const F a) { return std::round(a); }
template <typename F> F ceil(const F a) { return std::ceil(a); }
template <typename F> F floor(const F a) { return std::floor(a); }
template <typename F> F sqrt(const F a) { return std::sqrt(a); }
template <typename F> F sin(const F a) { return std::sin(a); }
template <typename F> F cos(const F a) { return std::cos(a); }
template <typename F> F tan(const F a) { return std::tan(a); }
template <typename F> F atan(const F a) { return std::atan(a); }
template <typename F> F atan2(const F a, const F b) { return std::atan2(a, b); }
template <typename F> F pow(const F a, const F b) { return std::pow(a, b); }
template <typename F> F exp(const F a) { return std::exp(a); }
template <typename F> F log(const F a) { return std::log(a); }
template <typename F> constexpr F numeric_max() {
  return std::numeric_limits<F>::max();
}

#endif
}
