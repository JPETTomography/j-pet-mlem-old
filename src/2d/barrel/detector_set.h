#pragma once

#include <iostream>

#include "square_detector.h"
#include "circle_detector.h"
#include "util/array.h"
#include "lor.h"

#include "symmetry_descriptor.h"

#if !__CUDACC__
#include "util/svg_ostream.h"
#include <vector>  // multi-ring detector construction
#endif
#ifndef MAX_DETECTORS
#define MAX_DETECTORS 256
#endif

/// Two-dimensional PET
namespace PET2D {
/// Two-dimensional PET barrel
namespace Barrel {

template <typename D> class ScannerBuilder;

template <typename DetectorType, std::size_t MaxDet, typename SType>
class DetectorSet : public util::array<MaxDet, DetectorType> {
 public:
  using Detector = DetectorType;
  using S = SType;
  using F = typename Detector::F;
  using LOR = Barrel::LOR<S>;
  using Point = PET2D::Point<F>;
  using Vector = PET2D::Vector<F>;
  using Circle = PET2D::Circle<F>;
  using Base = util::array<MaxDet, Detector>;
  using CircleDetector = Barrel::CircleDetector<F>;
  using Indices = util::array<MaxDet, S>;

  static const size_t MaxDetectors = MaxDet;
  struct Response {
    LOR lor;
    F dl;
  };

  /// Makes an empty detector set.
  DetectorSet(F radius = 1, F outer_radius = F(1.5))
      : Base(),
        fov_radius_(radius / M_SQRT2),
        c_inner(radius),
        c_outer(outer_radius),
        n_symmetries_(0) {}

  enum class TestCase {
    TEST_8_SQUARE_DETECTORS,
  };

  /// Makes new detector using hardcoded test case
  DetectorSet(TestCase test_case,  ///< test case
              F radius,            ///< radius of ring
              F w_detector,        ///< width of single detector (along ring)
              F h_detector,        ///< height/depth of single detector
                                   ///< (perpendicular to ring)
              F d_detector = 0     ///< diameter of circle single detector is
                                   ///< inscribed in
              )
      : Base(),
        fov_radius_(radius / M_SQRT2),
        c_inner(radius),
        c_outer(radius + (d_detector > 0 ? d_detector : h_detector)) {

    Detector detector_base(w_detector, h_detector, d_detector);

    switch (test_case) {
      case TestCase::TEST_8_SQUARE_DETECTORS:
        this->push_back(detector_base + Point(-radius, -radius));
        this->push_back(detector_base + Point(-radius, 0));
        this->push_back(detector_base + Point(-radius, radius));
        this->push_back(detector_base + Point(radius, -radius));
        this->push_back(detector_base + Point(radius, 0));
        this->push_back(detector_base + Point(radius, radius));
        this->push_back(detector_base + Point(0, -radius));
        this->push_back(detector_base + Point(0, radius));
        break;
      default:
        throw("unknown test case");
    }
  }

  int n_symmetries() const { return n_symmetries_; }

  void set_fov_radius(F fov_radius) { this->fov_radius_ = fov_radius; }

  F radius() const { return c_inner.radius; }
  F outer_radius() const { return c_outer.radius; }
  F max_dl(F max_bias_size) const { return 2 * c_outer.radius + max_bias_size; }

  /// Quantizes position across lor
  _ static S quantize_tof_position(F position,    ///< position across lor
                                   F step_size,   ///< step size
                                   S n_positions  ///< number of positions
                                   ) {
    // number of positions if always even, lower half are negative positions
    // where 0 means position closests to detector with higher index
    // maximum means position closests to detector with lower index
    if (position < 0)
      return n_positions / 2 - 1 -
             static_cast<S>(compat::floor(-position / step_size));
    else
      return static_cast<S>(compat::floor(position / step_size)) +
             n_positions / 2;
  }

  /// Returns number of position steps (indexes)
  S n_tof_positions(F step_size,     ///< step size
                    F max_bias_size  ///< possible bias (fuzz) maximum size
                    ) const {
    // since position needs to be symmetric against (0,0) number must be even
    return (static_cast<S>(ceil(2 * max_dl(max_bias_size) / step_size)) + 1) /
           2 * 2;
  }

  SymmetryDescriptor<S>* symmetry_descriptor() const {
    return symmetry_descriptor_;
  }

#if !__CUDACC__
  friend util::svg_ostream<F>& operator<<(util::svg_ostream<F>& svg,
                                          DetectorSet& cd) {
    svg << cd.c_outer;
    svg << cd.c_inner;

    svg << "<g id=\"photomultipiers\">" << std::endl;
    for (auto& detector : cd.c_detectors) {
      svg << detector;
    }
    svg << "</g>" << std::endl;

    svg << "<g id=\"scintillators\">" << std::endl;

    for (auto& detector : cd) {
      svg << detector;
    }

    svg << "</g>" << std::endl;

    return svg;
  }

  void to_mathematica(std::ostream& m_out) const {
    int i = 0;
    std::string delimiter;
    m_out << "{\"Detector\"->";
    m_out << "{\n";

    for (auto& detector : (*this)) {
      m_out << delimiter << "{ " << i << " , ";
      detector.to_mathematica(m_out);
      m_out << "}\n";
      i++;
      delimiter = ",";
    }
    m_out << "},";
    m_out << "\"Symmetries\"->{";
    symmetry_descriptor()->to_mathematica(m_out);
    m_out << "}";
    m_out << "}";
  }

#endif

  void push_back(const Detector& detector) {
    Base::push_back(detector);
    c_detectors.push_back(this->back().circumscribe_circle());
  }

  void push_back(Detector&& detector) {
    Base::push_back(std::forward<Detector>(detector));
    c_detectors.push_back(this->back().circumscribe_circle());
  }

  template <class... Args> void emplace_back(Args&&... args) {
    Base::emplace_back(std::forward<Args>(args)...);
    c_detectors.push_back(this->back().circumscribe_circle());
  }

  const CircleDetector& circumscribed(int i) const { return c_detectors[i]; }

  _ F fov_radius() const { return fov_radius_; }

  template <typename D> friend class ScannerBuilder;

 protected:
  F fov_radius_;
  util::array<MaxDet, CircleDetector> c_detectors;
  Circle c_inner;
  Circle c_outer;

  int n_symmetries_;
  SymmetryDescriptor<S>* symmetry_descriptor_;

 private:
};
}
}
