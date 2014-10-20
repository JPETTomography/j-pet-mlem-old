#pragma once

#include <cstddef>
#include <vector>

#include "vector.h"

namespace Math {

// FIXME: we asssume that this is diagonal matrix for now
template <size_t Dim, typename FType = double>
class Matrix : std::vector<FType> {
 public:
  typedef FType F;
  typedef Math::Vector<Dim> Vector;

  Matrix(std::initializer_list<F> l) : std::vector<F>(l) {}

  F& operator()(size_t i, size_t j) { return (*this)[Dim * i + j]; }

  Vector operator*(const Vector& v) const {
    Vector r{};
    for (int i = 0; i < Dim; ++i) {
      r[i] = v[i] * (*this)(i, i);
    }
    return r;
  }
};
}  // Math
