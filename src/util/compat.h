#pragma once

namespace compat {

#if __CUDACC__

#define GPU 1
#define $ __device__ __host__

template <typename F> F min(const F a, const F b);
template <> float min<float>(const float a, const float b) $ {
  return fminf(a, b);
}

template <typename F> F max(const F a, const F b);
template <> float max<float>(const float a, const float b) $ {
  return fmaxf(a, b);
}

template <typename F> F sqrt(const F a);
template <> float sqrt<float>(const float a) $ { return sqrtf(a); }

template <typename F> F sin(const F a);
template <> float sin<float>(const float a) $ { return sinf(a); }

template <typename F> F cos(const F a);
template <> float cos<float>(const float a) $ { return cosf(a); }

template <typename F> F tan(const F a);
template <> float tan<float>(const float a) $ { return tanf(a); }

template <typename F> F atan(const F a);
template <> float atan<float>(const float a) $ { return atanf(a); }

#else

#define $

template <typename F> F min(const F a, const F b) { return std::min(a, b); }
template <typename F> F max(const F a, const F b) { return std::max(a, b); }
template <typename F> F sqrt(const F a) { return std::sqrt(a); }
template <typename F> F sin(const F a) { return std::sin(a); }
template <typename F> F cos(const F a) { return std::cos(a); }
template <typename F> F tan(const F a) { return std::tan(a); }
template <typename F> F atan(const F a) { return std::atan(a); }

#endif
}
