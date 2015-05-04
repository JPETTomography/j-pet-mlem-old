#pragma once
#include <iostream>

#include <initializer_list>

#include "3d/geometry/vector.h"

namespace PET3D {

template <typename FType> class Matrix {
 public:
  using F = FType;

  static Matrix identity() {
    Matrix mat;
    mat(0, 0) = 1;
    mat(1, 1) = 1;
    mat(2, 2) = 1;
    return mat;
  }

  Matrix() : rep_() {}

  Matrix(std::initializer_list<F> elements) {
    int i = 0;
    for (auto it = elements.begin(); it != elements.end(); it++, i++) {
      rep_[i] = *it;
    }
  }

  F& operator()(int i) { return rep_[i]; }
  F operator()(int i) const { return rep_[i]; }

  F& operator()(int i, int j) { return rep_[i * 3 + j]; }
  F operator()(int i, int j) const { return rep_[i * 3 + j]; }

 private:
  F rep_[3 * 3];
};

template <typename FType>
Vector<FType> operator*(Matrix<FType> mat, Vector<FType> vec) {
  Vector<FType> result;

  result.x = mat(0, 0) * vec.x + mat(0, 1) * vec.y + mat(0, 2) * vec.z;
  result.y = mat(1, 0) * vec.x + mat(1, 1) * vec.y + mat(1, 2) * vec.z;
  result.z = mat(2, 0) * vec.x + mat(2, 1) * vec.y + mat(2, 2) * vec.z;

  return result;
}
}
