#pragma once

#include <ostream>
#include <random>

#include "2d/barrel/generic_scanner.h"
#include "3d/geometry/event.h"
#include "cmdline.h"
#include "3d/hybrid/options.h"

/// Three-dimensional PET
namespace PET3D {
/// Three-dimensional PET with longitudinal direction(z) added to the barrel
/// detector
namespace Hybrid {

/// 3D scanner made of several scintillators

template <typename Scanner2D> class Scanner {
 public:
  using F = typename Scanner2D::F;
  using S = typename Scanner2D::S;
  using Point = PET3D::Point<F>;
  using Point2D = PET2D::Point<F>;
  using Vector = PET3D::Vector<F>;
  using Vector2D = PET2D::Vector<F>;
  using Event = PET3D::Event<F>;
  using LOR = typename Scanner2D::LOR;
  using Indices = typename Scanner2D::Indices;
  using Event2D = typename Scanner2D::Event;

  static Scanner build_scanner_from_cl(const cmdline::parser& cl) {
    try {
      Scanner2D barrel =
          PET2D::Barrel::ScannerBuilder<Scanner2D>::build_multiple_rings(
              PET3D_LONGITUDINAL_SCANNER_CL(cl, F));
      barrel.set_fov_radius(cl.get<double>("fov-radius"));
      return Scanner(barrel, F(cl.get<double>("length")));
    } catch (const char* msg) {
      std::cerr << msg << "\n";
      exit(-1);
    }
  }

  /// Scanner full response
  ///
  /// Contains the full(redundant) information information about event and
  /// scanner response.
  struct FullResponse {
    S detector1, detector2;
    Point d1_entry, d1_exit, d1_deposition;
    Point d2_entry, d2_exit, d2_deposition;
    Point origin;

    friend std::ostream& operator<<(std::ostream& out,
                                    const FullResponse& response) {

      out << response.detector1 << " " << response.detector2;
      out << " " << response.d1_entry << " " << response.d1_exit << " "
          << response.d1_deposition;
      out << " " << response.d2_entry << " " << response.d2_exit << " "
          << response.d2_deposition;
      return out;
    };
  };
  /// Scanner response
  ///
  /// Represents information actually detected by the scanner on single event.
  struct Response {
    LOR lor;
    F z_up;
    F z_dn;
    F dl;

    friend std::ostream& operator<<(std::ostream& out,
                                    const Response& response) {
      out << (int)response.lor.first << " " << (int)response.lor.second << " ";
      out << response.z_up << " " << response.z_dn << " " << response.dl;
      return out;
    }
  };

  Scanner(const Scanner2D& barrel, F length)
      : barrel(barrel),
        length(length),
        half_length(length / 2),
        sigma_z_(),
        sigma_dl_() {}

  void set_sigmas(F sigma_z, F sigma_dl) {
    sigma_z_ = sigma_z;
    sigma_dl_ = sigma_dl;
  }

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
                       FullResponse& response   ///< response (LOR+zu+zd+dl)
                       ) const {

    if (escapes_through_endcap(e))
      return 0;

    Indices left, right;
    Event2D event_xy = e.to_barrel_event();

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

    if (detector1 > detector2) {
      response.detector1 = detector1;
      response.detector2 = detector2;

      response.d1_entry = d1_p1;
      response.d1_exit = d1_p2;
      response.d1_deposition = d1_deposit;

      response.d2_entry = d2_p1;
      response.d2_exit = d2_p2;
      response.d2_deposition = d2_deposit;
    } else {
      response.detector1 = detector2;
      response.detector2 = detector1;

      response.d1_entry = d2_p1;
      response.d1_exit = d2_p2;
      response.d1_deposition = d2_deposit;

      response.d2_entry = d1_p1;
      response.d2_exit = d1_p2;
      response.d2_deposition = d1_deposit;
    }
    response.origin = e.origin;

    return 2;
  }

  Response noErrorResponse(const FullResponse& full_response) const {
    Response response;
    F length1 = (full_response.d1_deposition - full_response.origin).length();
    F length2 = (full_response.d2_deposition - full_response.origin).length();

    response.lor = LOR(full_response.detector1, full_response.detector2);

    response.dl = length1 - length2;
    response.z_up = full_response.d1_deposition.z;
    response.z_dn = full_response.d2_deposition.z;

    return response;
  }

  template <typename RNG>
  Response errorResponse(RNG& rng, const FullResponse& full_response) const {
    std::normal_distribution<F> dist_z(0, sigma_z_);
    std::normal_distribution<F> dist_dl(0, sigma_dl_);
    Response response = noErrorResponse(full_response);

    response.z_up += dist_z(rng);
    response.z_dn += dist_z(rng);
    response.dl += dist_dl(rng);

    return response;
  }

  template <class RandomGenerator, class AcceptanceModel>
  _ short detect(RandomGenerator& gen,    ///< random number generator
                 AcceptanceModel& model,  ///< acceptance model
                 const Event& e,          ///< event to be detected
                 Response& response       ///< response (LOR+zu+zd+dl)
                 ) const {

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
                        Event2D e_xy,
                        S& detector,
                        Point& entry,
                        Point& exit,
                        Point& deposit) const {

    for (auto i : indices) {
      Point2D p1_xy, p2_xy;
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

  F sigma_z() const { return sigma_z_; }
  F sigma_dl() const { return sigma_dl_; }

  const Scanner2D barrel;
  const F length;
  const F half_length;

 private:
  F sigma_z_;
  F sigma_dl_;
};

}  // Longitudinal

}  // PET3D
