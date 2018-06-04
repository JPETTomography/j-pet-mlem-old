#pragma once

#include <type_traits>

#include "util/cuda/compat.h"
#include "util/read.h"

namespace PET3D {

/// Discreete coordinates 3D voxel
template <typename SType> struct Voxel {
  using S = SType;
  using Size = typename std::common_type<S, int>::type;

  _ Voxel(S x, S y, S z) : x(x), y(y), z(z) {}
  _ Voxel() = default;

  S x, y, z;

#if !__CUDACC__
  /// Constructs Voxel from stream
  Voxel(std::istream& in)
      : x(util::read<S>(in)), y(util::read<S>(in)), z(util::read<S>(in)) {}
#endif

  /// Index for given width & depth
  _ Size index(S width, S height) const {
    return (static_cast<Size>(z) * height + y) * width + x;
  }

  _ bool operator!=(const Voxel& v) const { return x != v.x || y != v.y; }

  _ bool operator==(const Voxel& v) const { return x == v.x && y == v.y; }

  _ bool operator<(const Voxel& v) const {
    return z < v.z || (z == v.z && (y < v.y || (y == v.y && x < v.x)));
  }

  _ void clamp(const Voxel& tl, const Voxel& br) {
    x = compat::min(br.x, compat::max(tl.x, x));
    y = compat::min(br.y, compat::max(tl.y, y));
    z = compat::min(br.z, compat::max(tl.z, z));
  }

  _ Size distance_from_origin2() const {
    return static_cast<Size>(x) * x + static_cast<Size>(y) * y * z;
  }

  _ Voxel& operator+=(const Voxel& rhs) {
    x += rhs.x, y += rhs.y, z += rhs.z;
    return *this;
  }

  _ Voxel& operator-=(const Voxel& rhs) {
    x -= rhs.x, y -= rhs.y, z -= rhs.z;
    return *this;
  }

  _ Voxel operator+(const Voxel& rhs) const { return Voxel(*this) += rhs; }

  _ Voxel operator-(const Voxel& rhs) const { return Voxel(*this) -= rhs; }
};

}  // PET3D

#ifdef TEST_CASE
namespace Catch {
template <typename SType> struct StringMaker<PET3D::Voxel<SType>> {
  static std::string convert(const PET3D::Voxel<SType>& v) {
    std::ostringstream oss;
    oss << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    return oss.str();
  }
};
}
#endif
