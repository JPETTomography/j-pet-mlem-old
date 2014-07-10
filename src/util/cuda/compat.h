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

namespace compat {

#if __CUDACC__

#define _ __device__ __host__
#define constexpr

  _ float min(const float a, const float b) { return fminf(a, b); }
  _ float max(const float a, const float b) { return fmaxf(a, b); }
  _ float ceil(const float a) { return ceilf(a); }
  _ float floor(const float a) { return floorf(a); }
  _ float sqrt(const float a) { return sqrtf(a); }
  _ float sin(const float a) { return sinf(a); }
  _ float cos(const float a) { return cosf(a); }
  _ float tan(const float a) { return tanf(a); }
  _ float atan(const float a) { return atanf(a); }
  _ float pow(const float a, const float b) { return powf(a, b); }
  _ float exp(const float a) { return expf(a); }

#else

#define _

#if _MSC_VER
  template <typename F> F min(const F a, const F b) { return (a < b) ? a : b; }
  template <typename F> F max(const F a, const F b) { return (a > b) ? a : b; }
#else
  template <typename F> F min(const F a, const F b) { return std::min(a, b); }
  template <typename F> F max(const F a, const F b) { return std::max(a, b); }
#endif
  template <typename F> F ceil(const F a) { return std::ceil(a); }
  template <typename F> F floor(const F a) { return std::floor(a); }
  template <typename F> F sqrt(const F a) { return std::sqrt(a); }
  template <typename F> F sin(const F a) { return std::sin(a); }
  template <typename F> F cos(const F a) { return std::cos(a); }
  template <typename F> F tan(const F a) { return std::tan(a); }
  template <typename F> F atan(const F a) { return std::atan(a); }
  template <typename F> F pow(const F a, const F b) { return std::pow(a, b); }
  template <typename F> F exp(const F a) { return std::exp(a); }

#endif
}