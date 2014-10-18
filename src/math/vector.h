#pragma once

#include <cstddef>
#include <vector>

namespace Math {

template <size_t Dim, typename FType = double>
class Vector : std::vector<FType> {
 public:
  typedef FType F;

  Vector(std::initializer_list<F> l) : std::vector<F>(l) {}

  Vector& operator*=(const Vector& v) {
    for (size_t i = 0; i < Dim; ++i) {
      (*this)[i] *= v[i];
    }
    return *this;
  }

  Vector operator*(const Vector& v) const {
    Vector r(*this);
    r *= v;
    return r;
  }

  Vector operator*(const Vector&& v) const {
    v *= *this;
    return v;
  }
};
}  // Math
