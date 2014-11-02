#pragma once

#include "catch.hpp"

inline Approx operator"" _e13(long double value) {
  return Approx(value).epsilon(1e-13);
}

class InverseApprox : public Approx {
  using Approx::Approx;
  friend bool operator==(double lhs, InverseApprox const& rhs) {
    return operator==(-lhs, static_cast<Approx>(rhs));
  }
  friend bool operator==(InverseApprox const& lhs, double rhs) {
    return operator==(rhs, lhs);
  }
  friend bool operator!=(double lhs, InverseApprox const& rhs) {
    return !operator==(lhs, rhs);
  }
  friend bool operator!=(InverseApprox const& lhs, double rhs) {
    return !operator==(rhs, lhs);
  }
};
inline InverseApprox& operator-(Approx& other) {
  return static_cast<InverseApprox&>(other);
}
inline InverseApprox&& operator-(Approx&& other) {
  return static_cast<InverseApprox&&>(other);
}
