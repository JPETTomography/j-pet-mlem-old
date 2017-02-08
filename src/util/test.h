#pragma once

#include <cstdlib>

#include "catch.hpp"

inline Approx operator"" _e13(long double value) {
  return Approx(value).epsilon(1e-13);
}

inline Approx operator"" _e7(long double value) {
  return Approx(value).epsilon(1e-7);
}

#if !defined(DIR_SEP) && (defined(_WIN32) || defined(WIN32))
#define DIR_SEP '\\'
#else
#define DIR_SEP '/'
#endif

inline std::string operator"" _temp(const char* base_name, size_t len) {
  const char* tmpdir = std::getenv("TMPDIR");
  if (!tmpdir)
    tmpdir = std::getenv("TEMP");
  if (!tmpdir)
    tmpdir = std::getenv("TMP");
#if __linux
  if (tmpdir == nullptr)
    tmpdir = "/tmp";
#endif
  REQUIRE(tmpdir != nullptr);
  return std::string(tmpdir) + DIR_SEP + std::string(base_name, len);
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

template <typename Vec> class VecApprox {
 public:
  explicit VecApprox(const Vec& value)
      : m_epsilon(std::numeric_limits<float>::epsilon() * 100),
        m_scale(1.0),
        m_value(value) {}

  VecApprox(VecApprox const& other)
      : m_epsilon(other.m_epsilon),
        m_scale(other.m_scale),
        m_value(other.m_value) {}

  static VecApprox custom() { return VecApprox(Vec(0, 0, 0)); }

  VecApprox operator()(double value) {
    VecApprox approx(value);
    approx.epsilon(m_epsilon);
    approx.scale(m_scale);
    return approx;
  }

  friend bool operator==(Vec lhs, VecApprox const& rhs) {
    // Thanks to Richard Harris for his help refining this formula
    bool x_check =
        fabs(lhs.x - rhs.m_value.x) <
        rhs.m_epsilon *
            (rhs.m_scale + (std::max)(fabs(lhs.x), fabs(rhs.m_value.x)));
    bool y_check =
        fabs(lhs.y - rhs.m_value.y) <
        rhs.m_epsilon *
            (rhs.m_scale + (std::max)(fabs(lhs.y), fabs(rhs.m_value.y)));
    bool z_check =
        fabs(lhs.z - rhs.m_value.z) <
        rhs.m_epsilon *
            (rhs.m_scale + (std::max)(fabs(lhs.z), fabs(rhs.m_value.z)));

    return x_check && y_check && z_check;
  }

  friend bool operator==(VecApprox const& lhs, Vec rhs) {
    return operator==(rhs, lhs);
  }

  Approx& epsilon(double newEpsilon) {
    m_epsilon = newEpsilon;
    return *this;
  }

  Approx& scale(double newScale) {
    m_scale = newScale;
    return *this;
  }

  std::string toString() const {
    std::ostringstream oss;
    oss << "Approx( " << m_value << " )";
    return oss.str();
  }

  template <typename V> friend VecApprox<V> VApprox(V vec);

 private:
  double m_epsilon;
  double m_scale;
  Vec m_value;
};

template <typename V> VecApprox<V> VApprox(V vec) { return VecApprox<V>(vec); }

namespace Catch {
template <typename Vec> struct StringMaker<VecApprox<Vec>> {
  static std::string convert(const VecApprox<Vec>& vec) {
    return vec.toString();
  }
};
}

// Stupid but working workaround for making syntax highlight working for tests
// in QtCreator. Use TEST("...") instead TEST_CASE("...").
#define TEST TEST_CASE

/// \endcond
