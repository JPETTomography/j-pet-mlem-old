#pragma once
#include <iostream>

#include <initializer_list>

#include "3d/geometry/vector.h"

namespace PET3D {

/// 3D (3x3) matrix stored in row major manner
template <typename FType> class Matrix {
 public:
  using F = FType;
  using Vector = PET3D::Vector<F>;

  static Matrix identity() {
    Matrix mat;
    mat(0, 0) = 1;
    mat(1, 1) = 1;
    mat(2, 2) = 1;
    return mat;
  }

  static Matrix skew(const Vector& v) {
    Matrix mat;
    mat(0, 1) = -v.z;
    mat(0, 2) = v.y;
    mat(1, 0) = v.z;
    mat(1, 2) = -v.x;
    mat(2, 0) = -v.y;
    mat(2, 1) = v.x;

    return mat;
  }

  static Matrix rotate(const Vector& u, F angle) {
    auto s = compat::sin(angle);
    auto c = compat::cos(angle);
    auto t = 1 - c;
    // clang-format off
    return Matrix{ u.x*u.x*t + c,     u.x*u.y*t - u.z*s,  u.x*u.z*t + u.y*s,
                   u.y*u.x*t + u.z*s, u.y*u.y*t + c,      u.y*u.z*t - u.x*s,
                   u.z*u.x*t - u.y*s, u.z*u.y*t + u.x*s,  u.z*u.z*t + c      };
    // clang-format on
  }

  Matrix() : rep_() {}

  Matrix(std::initializer_list<F> elements) {
    int i = 0;
    for (auto it = elements.begin(); it != elements.end(); it++, i++) {
      rep_[i] = *it;
    }
  }

  template <typename I> Matrix(I begin, I end) {
    int i = 0;
    for (auto it = begin; it != end; it++, i++) {
      rep_[i] = *it;
    }
  }

  F& operator()(int i) { return rep_[i]; }
  F operator()(int i) const { return rep_[i]; }

  F& operator()(int row, int col) { return rep_[row * 3 + col]; }
  F operator()(int row, int col) const { return rep_[row * 3 + col]; }

  Matrix& operator+=(const Matrix& rhs) {
    for (int i = 0; i < 9; i++)
      rep_[i] += rhs(i);
    return *this;
  }

  Matrix& operator-=(const Matrix& rhs) {
    for (int i = 0; i < 9; i++)
      rep_[i] -= rhs(i);
    return *this;
  }

  Matrix& operator*=(F s) {
    for (int i = 0; i < 9; i++)
      rep_[i] *= s;
    return *this;
  }

 private:
  F rep_[3 * 3];
};

template <typename FType>
Matrix<FType> operator+(Matrix<FType> lhs, Matrix<FType> rhs) {
  Matrix<FType> res(lhs);
  res += rhs;
  return res;
}

template <typename FType>
Matrix<FType> operator-(Matrix<FType> lhs, Matrix<FType> rhs) {
  Matrix<FType> res(lhs);
  res -= rhs;
  return res;
}

template <typename FType>
Vector<FType> operator*(Matrix<FType> mat, Vector<FType> vec) {
  Vector<FType> result;

  result.x = mat(0, 0) * vec.x + mat(0, 1) * vec.y + mat(0, 2) * vec.z;
  result.y = mat(1, 0) * vec.x + mat(1, 1) * vec.y + mat(1, 2) * vec.z;
  result.z = mat(2, 0) * vec.x + mat(2, 1) * vec.y + mat(2, 2) * vec.z;

  return result;
}

template <typename FType> Matrix<FType> operator*(Matrix<FType> mat, FType s) {
  Matrix<FType> result(mat);
  result *= s;
  return result;
}

template <typename FType> Matrix<FType> operator*(FType s, Matrix<FType> mat) {
  Matrix<FType> result(mat);
  result *= s;
  return result;
}

template <typename FType> Matrix<FType> transpose(const Matrix<FType>& mat) {
  Matrix<FType> result;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      result(i, j) = mat(j, i);

  return result;
}

}  // PET3D
