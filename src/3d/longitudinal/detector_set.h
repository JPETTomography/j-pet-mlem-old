#include<iostream>

#include "2d/barrel/detector_set.h"
#include "3d/geometry/event.h"

/// Three-dimensional PET
namespace PET3D {
/// Three-dimensional PET with longitudinal direction(z) added to the barrel
/// detector
namespace Longitudinal {

/// 3D detector made of several scintillators

template <typename DetectorSet2D,
          std::size_t MaxDetectors,
          typename FType = float,
          typename SType = int>
class DetectorSet {
 public:
  using Event = PET3D::Event<FType>;
  using BarrelEvent = PET2D::Barrel::Event<FType>;

  DetectorSet(const DetectorSet2D& barrel_a, FType length_a)
      : barrel(barrel_a), length(length_a), half_length(length / 2) {}

  bool escapes_through_endcap(const Event& event) {
    if (event.direction.z == 0.0)
      return false;

    FType r2 = barrel.radius() * barrel.radius();

    {
      FType t_right = (half_length - event.origin.z) / event.direction.z;
      FType x_right = event.origin.x + t_right * event.direction.x;
      FType y_right = event.origin.y + t_right * event.direction.y;
      if (x_right * x_right + y_right * y_right < r2)
        return true;
    }
    {
      FType t_left = -(half_length + event.origin.z) / event.direction.z;
      FType x_left = event.origin.x + t_left * event.direction.x;
      FType y_left = event.origin.y + t_left * event.direction.y;
      if (x_left * x_left + y_left * y_left < r2)
        return true;
    }

    return false;
  }

 private:
  const DetectorSet2D& barrel;
  const FType length;
  const FType half_length;
};

}  // Longitudinal
}  // PET3D
