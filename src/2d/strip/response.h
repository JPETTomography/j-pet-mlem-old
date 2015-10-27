#pragma once

#if !__CUDACC__
#include <ostream>
#endif

#include "util/cuda/compat.h"

namespace PET2D {
namespace Strip {

/// Strip scanner response
template <typename FType> struct Response {
  using F = FType;
  F z_u;
  F z_d;
  F dl;

  _ Response(F z_u, F z_d, F dl) : z_u(z_u), z_d(z_d), dl(dl) {}
  _ Response() {}

  _ void transform(F R, F& tan, F& y, F& z) const {
    tan = this->tan(R);
    y = this->y(tan);
    z = this->z(y, tan);
  }

 private:
  _ F tan(const F R) const { return (z_u - z_d) / (2 * R); }

  _ F y(const F tan) const {
    return -F(0.5) * (dl / compat::sqrt(1 + tan * tan));
  }

  _ F z(const F y, const F tan) const {
    return F(0.5) * (z_u + z_d + (2 * y * tan));
  }

#if !__CUDACC__
  friend std::ostream& operator<<(std::ostream& out, const Response& response) {
    out << response.z_u << ' ' << response.z_d << ' ' << response.dl;
    return out;
  }
#endif
};

}  // Strip
}  // PET2D
