#pragma once

#include <iostream>

#include "2d/barrel/detector_set.h"
#include "3d/geometry/event.h"

/// Three-dimensional PET
namespace PET3D {
/// Three-dimensional PET with longitudinal direction(z) added to the barrel
/// detector
namespace Longitudinal {

/// 3D detector made of several scintillators

template <typename DetectorSet2D> class DetectorSet {
 public:
  using F = typename DetectorSet2D::F;
  using S = typename DetectorSet2D::S;
  using Point = PET3D::Point<F>;
  using Point2D = PET2D::Point<F>;
  using Vector = PET3D::Vector<F>;
  using Vector2D = PET2D::Vector<F>;
  using Event = PET3D::Event<F>;
  using LOR = typename DetectorSet2D::LOR;
  using Indices = typename DetectorSet2D::Indices;
  using BarrelEvent = typename DetectorSet2D::Event;
  using BarrelType = DetectorSet2D;

  /**
 * @brief The FullResponse struct
 *
 * contains the full(redundant) information information about event and detector
 * response.
 */
  struct FullResponse {
    S detector1, detector2;
    Point d1_entry, d1_exit, d1_deposition;
    Point d2_entry, d2_exit, d2_deposition;
    Point origin;
  };

  struct Response {
    LOR lor;
    F z_up;
    F z_dn;
    F dl;
  };

  DetectorSet(const DetectorSet2D& barrel, F length)
      : barrel(barrel), length(length), half_length(length / 2) {}

  bool escapes_through_endcap(const Event& event) const {
    if (event.direction.z == 0.0)
      return false;

    F r2 = barrel.radius() * barrel.radius();

    F t_right = (half_length - event.origin.z) / event.direction.z;
    F x_right = event.origin.x + t_right * event.direction.x;
    F y_right = event.origin.y + t_right * event.direction.y;
    if (x_right * x_right + y_right * y_right < r2)
      return true;

    F t_left = -(half_length + event.origin.z) / event.direction.z;
    F x_left = event.origin.x + t_left * event.direction.x;
    F y_left = event.origin.y + t_left * event.direction.y;
    if (x_left * x_left + y_left * y_left < r2)
      return true;

    return false;
  }

  template <class RandomGenerator, class AcceptanceModel>
  _ short exact_detect(RandomGenerator& gen,    ///< random number generator
                       AcceptanceModel& model,  ///< acceptance model
                       const Event& e,          ///< event to be detected
                       FullResponse& response) const {

    if (escapes_through_endcap(e))
      return 0;

    Indices left, right;
    BarrelEvent event_xy = e.to_barrel_event();

    barrel.close_indices(event_xy, left, right);

    S detector1, detector2;

    Point d1_p1, d1_p2, d2_p1, d2_p2;
    Point d1_deposit, d2_deposit;

    if (!check_for_hits(gen,
                        model,
                        left,
                        e,
                        -1.0,
                        event_xy,
                        detector1,
                        d1_p1,
                        d1_p2,
                        d1_deposit) ||
        !check_for_hits(gen,
                        model,
                        right,
                        e,
                        1.0f,
                        event_xy,
                        detector2,
                        d2_p1,
                        d2_p2,
                        d2_deposit))
      return 0;

    response.detector1 = detector1;
    response.detector2 = detector2;

    response.d1_entry = d1_p1;
    response.d1_exit = d1_p2;
    response.d1_deposition = d1_deposit;

    response.d2_entry = d2_p1;
    response.d2_exit = d2_p2;
    response.d2_deposition = d2_deposit;
    response.origin = e.origin;

    return 2;
  }

  Response noErrorResponse(const FullResponse& full_response) const {
    Response response;
    F length1 = (full_response.d1_deposition - full_response.origin).length();
    F length2 = (full_response.d2_deposition - full_response.origin).length();

    response.lor = LOR(full_response.detector1, full_response.detector2);

    if (full_response.detector1 > full_response.detector2) {
      response.dl = length1 - length2;
      response.z_up = full_response.d1_deposition.z;
      response.z_dn = full_response.d2_deposition.z;
    } else {
      response.dl = length2 - length1;
      response.z_up = full_response.d2_deposition.z;
      response.z_dn = full_response.d1_deposition.z;
    }
    return response;
  }

  template <class RandomGenerator, class AcceptanceModel>
  _ short detect(RandomGenerator& gen,    ///< random number generator
                 AcceptanceModel& model,  ///< acceptance model
                 const Event& e,          ///< event to be detected
                 Response& response) const {

    FullResponse full_response;
    if (exact_detect(gen, model, e, full_response) == 2) {

      response = noErrorResponse(full_response);
      return 2;
    }
    return 0;
  }

  void reconstruct_3d_intersection_points(const Event& event,
                                          F dir,
                                          const Point2D& entry_xy,
                                          const Point2D& exit_xy,
                                          Point& entry,
                                          Point& exit) const {
    Vector2D dir_xy = event.direction.xy();
    Point2D origin_xy = event.origin.xy();
    F dz_over_dx_dxy = dir * event.direction.z / dir_xy.length();

    Vector2D displacement_entry = entry_xy - origin_xy;
    F displacement_entry_length = displacement_entry.length();
    F dz_entry = dz_over_dx_dxy * displacement_entry_length;

    entry.x = entry_xy.x;
    entry.y = entry_xy.y;
    entry.z = dz_entry + event.origin.z;

    Vector2D displacement_exit = exit_xy - origin_xy;
    F displacement_exit_length = displacement_exit.length();
    F dz_exit = dz_over_dx_dxy * displacement_exit_length;

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
                     F& depth) const {
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
                        F dir,
                        BarrelEvent e_xy,
                        S& detector,
                        Point& entry,
                        Point& exit,
                        Point& deposit) const {

    for (auto i : indices) {
      Point2D p1_xy, p2_xy;
      Point2D origin_2d = e_xy;
      F depth;
      if (barrel.did_intersect(e_xy, i, p1_xy, p2_xy)) {
        reconstruct_3d_intersection_points(e, dir, p1_xy, p2_xy, entry, exit);
        if (entry.z > half_length || entry.z < -half_length)
          return false;
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

  const DetectorSet2D barrel;
  const F length;
  const F half_length;
};

}  // Longitudinal
}  // PET3D
