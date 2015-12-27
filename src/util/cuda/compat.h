#pragma once

#if !__CUDACC__
#include <cmath>
#include <algorithm>
#include <limits>
#else
#include <limits.h>
#include <float.h>
#endif

/// CUDA generic code compatibility macros & functions
////
/// Author
/// ------
/// Adam Strzelecki <adam.strzelecki@uj.edu.pl>
///
/// Discussion
/// ----------
/// Purpose of this header is to provide maximum share-ability between generic
/// CPU implementation and CUDA implementation. This is done using following:
///
/// 1. All mathematical function are exposed via compat namespace i.e.
///    compat::cos which map to stdlib or CUDA depending on build context.
///
/// 2. Special underscore macro `_` is used to mark functions and methods
///    compatible with CUDA and it is replaced with `__device__ __host__` when
///    building using CUDA.

namespace compat {

#if __CUDACC__

#define _ __device__ __host__

_ static inline int abs(const int a) { return ::abs(a); }
_ static inline int min(const int a, const int b) { return ::min(a, b); }
_ static inline int max(const int a, const int b) { return ::max(a, b); }
_ static inline float min(const float a, const float b) { return fminf(a, b); }
_ static inline float max(const float a, const float b) { return fmaxf(a, b); }
_ static inline float round(const float a) { return roundf(a); }
_ static inline float ceil(const float a) { return ceilf(a); }
_ static inline float floor(const float a) { return floorf(a); }
_ static inline float sqrt(const float a) { return sqrtf(a); }
_ static inline float sin(const float a) { return sinf(a); }
_ static inline float cos(const float a) { return cosf(a); }
_ static inline float tan(const float a) { return tanf(a); }
_ static inline float atan(const float a) { return atanf(a); }
_ static inline float atan2(const float a, const float b) {
  return atan2f(a, b);
}
_ static inline float pow(const float a, const float b) { return powf(a, b); }
_ static inline float exp(const float a) { return expf(a); }
_ static inline float log(const float a) { return logf(a); }

template <typename F> _ constexpr F numeric_max();
template <> _ constexpr float numeric_max<float>() { return FLT_MAX; }
template <> _ constexpr int numeric_max<int>() { return INT_MAX; }
template <> _ constexpr unsigned int numeric_max<unsigned int>() {
  return UINT_MAX;
}
template <> _ constexpr unsigned long numeric_max<unsigned long>() {
  return ULONG_MAX;
}

#else  // __CUDACC__

#define _

#if _MSC_VER
template <typename F> F abs(const F a) { return (a >= 0) ? a : -a; }
template <typename F> F min(const F a, const F b) { return (a < b) ? a : b; }
template <typename F> F max(const F a, const F b) { return (a > b) ? a : b; }
#else
template <typename F> F abs(const F a) { return std::abs(a); }
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

/// \cond PRIVATE

#endif  // __CUDACC__

/// \endcond

/// Drop-in replacement for `std::swap`
template <typename Assignable1, typename Assignable2>
_ inline void swap(Assignable1& a, Assignable2& b) {
  Assignable1 temp = a;
  a = b;
  b = temp;
}

}  // compat
