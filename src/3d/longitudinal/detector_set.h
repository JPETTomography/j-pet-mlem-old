#include <iostream>

#include "2d/barrel/detector_set.h"
#include "3d/geometry/event.h"

/// Three-dimensional PET
namespace PET3D {
/// Three-dimensional PET with longitudinal direction(z) added to the barrel
/// detector
namespace Longitudinal {

/// 3D detector made of several scintillators

template <typename DetectorSet2D,
        //  std::size_t MaxDetectors,
          typename FType = typename DetectorSet2D::F,
          typename SType = typename DetectorSet2D::S>
class DetectorSet {
 public:
  using Point = PET3D::Point<FType>;
  using Point2D = PET2D::Point<FType>;
  using Vector = PET3D::Vector<FType>;
  using Vector2D = PET2D::Vector<FType>;
  using Event = PET3D::Event<FType>;
  using LOR = typename DetectorSet2D::LOR;
  using Indices = typename DetectorSet2D::Indices;
  using BarrelEvent = typename DetectorSet2D::Event;

  struct Response {
    LOR lor;
    FType z_up;
    FType z_dn;
    FType dl;
  };

  DetectorSet(const DetectorSet2D& barrel, FType length)
      : barrel(barrel), length(length), half_length(length / 2) {}

  bool escapes_through_endcap(const Event& event) const {
    if (event.direction.z == 0.0)
      return false;

    FType r2 = barrel.radius() * barrel.radius();

    FType t_right = (half_length - event.origin.z) / event.direction.z;
    FType x_right = event.origin.x + t_right * event.direction.x;
    FType y_right = event.origin.y + t_right * event.direction.y;
    if (x_right * x_right + y_right * y_right < r2)
      return true;

    FType t_left = -(half_length + event.origin.z) / event.direction.z;
    FType x_left = event.origin.x + t_left * event.direction.x;
    FType y_left = event.origin.y + t_left * event.direction.y;
    if (x_left * x_left + y_left * y_left < r2)
      return true;

    return false;
  }

  template <class RandomGenerator, class AcceptanceModel>
  _ short detect(RandomGenerator& gen,    ///< random number generator
                 AcceptanceModel& model,  ///< acceptance model
                 const Event& e,          ///< event to be detected
                 Response& response) const {

    if (escapes_through_endcap(e))
      return 0;

    Indices left, right;
    BarrelEvent event_xy = e.to_barrel_event();

    barrel.close_indices(event_xy, left, right);

    SType detector1, detector2;

    Point d1_p1, d1_p2, d2_p1, d2_p2;
    Point deposit_d1, deposit_d2;

    if (!check_for_hits(gen,
                        model,
                        left,
                        e,
                        -1.0,
                        event_xy,
                        detector1,
                        d1_p1,
                        d1_p2,
                        deposit_d1) ||
        !check_for_hits(gen,
                        model,
                        right,
                        e,
                        1.0f,
                        event_xy,
                        detector2,
                        d2_p1,
                        d2_p2,
                        deposit_d2))
      return 0;

    response.lor = LOR(detector1, detector2);
    FType length1 = (deposit_d1 - e.origin).length();
    FType length2 = (deposit_d2 - e.origin).length();

    if (detector1 > detector2) {
      response.dl = length1 - length2;
      response.z_up = deposit_d1.z;
      response.z_dn = deposit_d2.z;
    } else {
      response.dl = length2 - length1;
      response.z_up = deposit_d2.z;
      response.z_dn = deposit_d1.z;
    }
    return 2;
  }

  void reconstruct_3d_intersection_points(const Event& event,
                                          FType dir,
                                          const Point2D& entry_xy,
                                          const Point2D& exit_xy,
                                          Point& entry,
                                          Point& exit) const {
    Vector2D dir_xy = event.direction.xy();
    Point2D origin_xy = event.origin.xy();
    FType dz_over_dx_dxy = dir * event.direction.z / dir_xy.length();

    Vector2D displacement_entry = entry_xy - origin_xy;
    FType displacement_entry_length = displacement_entry.length();
    FType dz_entry = dz_over_dx_dxy * displacement_entry_length;

    entry.x = entry_xy.x;
    entry.y = entry_xy.y;
    entry.z = dz_entry + event.origin.z;

    Vector2D displacement_exit = exit_xy - origin_xy;
    FType displacement_exit_length = displacement_exit.length();
    FType dz_exit = dz_over_dx_dxy * displacement_exit_length;

    exit.x = exit_xy.x;
    exit.y = exit_xy.y;
    exit.z = dz_exit + event.origin.z;

    if (displacement_entry_length > displacement_exit_length)
      std::swap(entry, exit);
  }

  template <class RandomGenerator, class AcceptanceModel>
  _ bool did_deposit(RandomGenerator& gen,
                     AcceptanceModel& model,
                     const Point& entry,
                     const Point& exit,
                     FType& depth) const {
    depth = model.deposition_depth(gen);
    if (depth * depth <= (exit - entry).length2())
      return true;
    return false;
  }

  template <class RandomGenerator, class AcceptanceModel>
  _ bool check_for_hits(RandomGenerator& gen,
                        AcceptanceModel& model,
                        const Indices& indices,
                        Event e,
                        FType dir,
                        BarrelEvent e_xy,
                        SType& detector,
                        Point& entry,
                        Point& exit,
                        Point& deposit) const {

    for (auto i : indices) {
      Point2D p1_xy, p2_xy;
      Point2D origin_2d = e_xy;
      FType depth;
      if (barrel.did_intersect(e_xy, i, p1_xy, p2_xy)) {
        reconstruct_3d_intersection_points(e, dir, p1_xy, p2_xy, entry, exit);
        if (did_deposit(gen, model, entry, exit, depth)) {
          Vector v = exit - entry;
          deposit = entry + v * (depth / v.length());
          detector = i;
          return true;
        }
      }
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
