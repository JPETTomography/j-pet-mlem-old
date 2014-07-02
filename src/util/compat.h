#pragma once

namespace compat {

#if __CUDACC__

#define $ __device__ __host__
#define constexpr

template <typename F> F min(const F a, const F b);
template <> float min<float>(const float a, const float b) $ {
  return fminf(a, b);
}

template <typename F> F max(const F a, const F b);
template <> float max<float>(const float a, const float b) $ {
  return fmaxf(a, b);
}

template <typename F> F ceil(const F a);
template <> float ceil<float>(const float a) $ { return ceilf(a); }

template <typename F> F floor(const F a);
template <> float floor<float>(const float a) $ { return floorf(a); }

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

template <typename F> F pow(const F a, const F b);
template <> float pow<float>(const float a, const float b) $ {
  return powf(a, b);
}

template <typename F> F exp(const F a);
template <> float exp<float>(const float a) $ { return expf(a); }

#else

#define $

template <typename F> F min(const F a, const F b) { return std::min(a, b); }
template <typename F> F max(const F a, const F b) { return std::max(a, b); }
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
