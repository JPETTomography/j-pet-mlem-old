#pragma once

#include <cstdlib>

#include "catch.hpp"

inline Approx operator"" _e13(long double value) {
  return Approx(value).epsilon(1e-13);
}

inline Approx operator"" _e7(long double value) {
  return Approx(value).epsilon(1e-7);
}

inline std::string operator"" _temp(const char* base_name, unsigned long len) {
  const char* tmpdir =
      std::getenv("TMPDIR") ?: std::getenv("TEMP") ?: std::getenv("TEMP");
#if __linux
  if (tmpdir == nullptr)
    tmpdir = "/tmp";
#endif
  REQUIRE(tmpdir != nullptr);
  return std::string(tmpdir) + "/" + std::string(base_name, len);
}

/// \cond PRIVATE

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

// Stupid but working workaround for making syntax highlight working for tests
// in QtCreator. Use TEST("...") instead TEST_CASE("...").
#define TEST TEST_CASE

/// \endcond
