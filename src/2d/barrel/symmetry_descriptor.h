#pragma once

#include "2d/geometry/pixel.h"
#include "2d/geometry/pixel_grid.h"

#if !__CUDACC__
#include "util/json.h"
#endif

namespace PET2D {
namespace Barrel {

enum Axis { X = 1, Y = 2, XY = 4 };

/// Describes detector ring symmetries and symmetry transformations
template <typename SType> class SymmetryDescriptor {
 public:
  using S = SType;

  constexpr static const S EIGHT = 8;

  SymmetryDescriptor(int n_detectors, int n_symmetries)
      : n_detectors(n_detectors), n_symmetries(n_symmetries) {
    detectors_ = new S[n_detectors * n_symmetries];
  }

  ~SymmetryDescriptor() { delete detectors_; }

  /// Symmetric detector
  S symmetric_detector(S detector, S symmetry) const {
    return detectors_[detector * EIGHT + symmetry];
  }

  // converts the pixels from the triangular cut
  // to full pixel_grid pixels
  template <typename FType>
  static Pixel<S> full_grid_pixel(const PET2D::PixelGrid<FType, S> grid,
                                  const Pixel<S>& pixel) {
    return Pixel<S>(pixel.x + grid.n_columns / 2, pixel.y + grid.n_rows / 2);
  }

  template <typename FType>
  static Pixel<S> symmetric_pixel(const PET2D::PixelGrid<FType, S> grid,
                                  const Pixel<S>& pixel,
                                  S symmetry) {
    Pixel<S> s_pixel(pixel);
    if (symmetry & Axis::X) {
      s_pixel.y = grid.n_row - s_pixel.y - 1;
    }
    if (symmetry & Axis::Y) {
      s_pixel.y = grid.n_columns - s_pixel.x - 1;
    }
    if (symmetry & Axis::XY) {
      std::swap(s_pixel.x, s_pixel.y);
    }
    return s_pixel;
  }

  /// Returns symmetric detector on a ring of n_detectors, assuming that
  /// detector zero is on the positive X-axis (rotation = 0).
  static S ring_symmetric_detector(S n_detectors, S detector, S symmetry) {
    if (symmetry & Axis::X) {
      detector = (n_detectors - detector) % n_detectors;
    }
    if (symmetry & Axis::Y) {
      detector = ((n_detectors + n_detectors / 2) - detector) % n_detectors;
    }
    if (symmetry & Axis::XY) {
      detector = ((n_detectors + n_detectors / 4) - detector) % n_detectors;
    }
    return detector;
  }

  /// Rotated symmetric detector
  ///
  /// Returns symmetric detector on a ring of n_detectors, assuming that
  /// detector zero is at the angle Pi/n_detectors (rotation = 0.5).
  static S rotated_ring_symmetric_detector(S n_detectors,
                                           S detector,
                                           S symmetry) {
    if (symmetry & Axis::X) {
      detector = (n_detectors - (detector + 1)) % n_detectors;
    }
    if (symmetry & Axis::Y) {
      detector =
          ((n_detectors + n_detectors / 2) - (detector + 1)) % n_detectors;
    }
    if (symmetry & Axis::XY) {
      detector =
          ((n_detectors + n_detectors / 4) - (detector + 1)) % n_detectors;
    }
    return detector;
  }

  void set_symmetric_detector(S detector, S symmetry, S symmetric_detector) {
    detectors_[detector * EIGHT + symmetry] = symmetric_detector;
  }

#if !__CUDACC__

  void serialize(std::ostream& out) {
    out << n_detectors << " " << n_symmetries << "\n";
    for (S d = 0; d < n_detectors; d++) {
      out << d << " ";
      for (S s = 0; s < n_symmetries; s++)
        out << detectors_[d * EIGHT + s] << " ";
      out << "\n";
    }
  }

  static SymmetryDescriptor* deserialize(std::istream& in) {
    S n_detectors, n_symmetries;
    in >> n_detectors >> n_symmetries;
    auto descriptor = new SymmetryDescriptor(n_detectors, n_symmetries);
    S dd;
    for (S d = 0; d < n_detectors; d++) {
      in >> dd;
      for (S s = 0; s < n_symmetries; s++) {
        S ds;
        in >> ds;
        descriptor->set_symmetric_detector(dd, s, ds);
      }
    }
    return descriptor;
  }

  operator json() const {
    json j;
    for (S i = 0; i < n_detectors; ++i) {
      json j_symmetric_detectors;
      for (S s = 0; s < n_symmetries; ++s) {
        j_symmetric_detectors.push_back(symmetric_detector(i, s));
      }
      j[i] = j_symmetric_detectors;
    }
    return j;
  }
#endif

  const S n_detectors;
  const S n_symmetries;

 private:
  S* detectors_;
};

}  // Barrel
}  // PET2D
