#pragma once

#include "square_detector.h"
#include "circle_detector.h"
#include "util/array.h"
#include "lor.h"

#include "symmetry_descriptor.h"

#if !__CUDACC__
#include "util/svg_ostream.h"
#include "util/json_ostream.h"
#endif

/// Two-dimensional PET
namespace PET2D {
/// Two-dimensional PET barrel
namespace Barrel {

template <class ScannerClass> class ScannerBuilder;

template <class DetectorClass,
          typename SType,
          std::size_t MaxDetetectorsSize = 192>
class DetectorSet : public util::array<MaxDetetectorsSize, DetectorClass> {
 public:
  using Detector = DetectorClass;
  using S = SType;
  static const size_t MaxDetectors = MaxDetetectorsSize;
  using F = typename Detector::F;
  using LOR = Barrel::LOR<S>;
  using Point = PET2D::Point<F>;
  using Vector = PET2D::Vector<F>;
  using Circle = PET2D::Circle<F>;
  using Base = util::array<MaxDetetectorsSize, Detector>;
  using CircleDetector = Barrel::CircleDetector<F>;
  using Indices = util::array<MaxDetetectorsSize, S>;
  using Event = Barrel::Event<F>;

  struct Response {
    LOR lor;
    F dl;
    int tof_position;

#if !__CUDACC__
    friend std::ostream& operator<<(std::ostream& out,
                                    const Response& response) {
      out << (int)response.lor.first << " " << (int)response.lor.second;
      out << " " << response.tof_position;
      out << " " << response.dl;
      return out;
    }
#endif
  };

  using FullResponse = Response;

  /// Makes an empty detector set.
  DetectorSet(F radius = 1, F outer_radius = F(1.5))
      : Base(),
        fov_radius_(radius / M_SQRT2),
        c_inner(radius),
        c_outer(outer_radius),
        n_symmetries_(0),
        tof_step_size_() {}

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

  F max_dl(F max_dl_error) const { return 2 * c_outer.radius + max_dl_error; }

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
  S n_tof_positions(F step_size,    ///< step size
                    F max_dl_error  ///< possible bias (fuzz) maximum size
                    ) const {
    // since position needs to be symmetric against (0,0) number must be even
    return (static_cast<S>(ceil(2 * max_dl(max_dl_error) / step_size)) + 1) /
           2 * 2;
  }

  SymmetryDescriptor<S>& symmetry_descriptor() const {
    return *symmetry_descriptor_;
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

  friend util::json_ostream& operator<<(util::json_ostream& json,
                                        const DetectorSet& ds) {
    bool next = false;

    json << "{\"Detector\":"
         << "[\n";

    for (auto& detector : ds) {
      json.delimiter(next) << detector;
    }

    json << "]\n"
         << "}";
    return json;
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

  template <class ScannerClass> friend class ScannerBuilder;

  _ void set_tof_step(F tof_step_size) { tof_step_size_ = tof_step_size; }
  _ F tof_step_size() const { return tof_step_size_; }

 protected:
  F fov_radius_;
  util::array<MaxDetetectorsSize, CircleDetector> c_detectors;
  Circle c_inner;
  Circle c_outer;
  int n_symmetries_;
  SymmetryDescriptor<S>* symmetry_descriptor_;

 private:
  F tof_step_size_;
};

}  // Barrel
}  // PET2D
