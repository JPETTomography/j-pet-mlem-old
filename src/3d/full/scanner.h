#ifndef SCANNER_H
#define SCANNER_H

#include <vector>

#include "3d/full/volume.h"
#include "3d/geometry/event.h"
#include "3d/geometry/ray.h"
#include "util/array.h"
#include "util/sort.h"

namespace PET3D {
namespace Full {

template <typename F, typename S, int MAX_VOLUMES = 2 << 9> class Scanner {
  using Event = PET3D::Event<F>;
  using Point = PET3D::Point<F>;
  using Ray = ray_tracing::Ray<F>;

  struct HalfResponse {
    S detector;
    Point entry, exit, deposition;
  };

 public:
  struct FullResponse {
    S detector1, detector2;
    Point d1_entry, d1_exit, d1_deposition;
    Point d2_entry, d2_exit, d2_deposition;
    Point origin;

    FullResponse() = default;
  };

  void add_volume(Volume<F>* vol) { volumes_.push_back(vol); };

  template <class RNG, class AcceptanceModel>
  S exact_detect(RNG& rng,                ///< random number generator
                 AcceptanceModel& model,  ///< acceptance model
                 const Event& event,      ///< event to be detected
                 FullResponse& response   ///< response (LOR+zu+zd+dl)
                 ) {

    response.origin = event.origin;

    Ray up(event.origin, event.direction);
    Ray dn(event.origin, -event.direction);
    S intersected_up = 0;
    using intersection_t = std::pair<int, ray_tracing::intersection_result<F>>;

    util::array<MAX_VOLUMES, intersection_t> intersected_volumes_up_;
    int hits = 0;
    int vol_up;
    for (int i = 0; i < volumes_.size(); i++) {
      auto v = volumes_[i];
      auto hit_up = v->intersects_with(up);
      if (hit_up.intersected) {

        intersected_volumes_up_.emplace_back(i, hit_up);
      }
    }

    util::heap_sort(intersected_volumes_up_.begin(),
                    intersected_volumes_up_.end(),
                    [&](const intersection_t& a, const intersection_t& b) {
                      return a.second.t_min < a.second.t_min;
                    });
    for (auto res : intersected_volumes_up_) {
      F l_up = res.second.t_max - res.second.t_min;
      F l_depth = model.deposition_depth(rng);
      if (l_depth < l_up) {
        hits++;
        response.detector1 = res.first;
        response.d1_entry = up(res.second.t_min);
        response.d1_deposition = up(res.second.t_min + l_depth);
        response.d1_exit = up(res.second.t_max);
        break;
      }
    }

    if (hits < 1)
      return hits;

    util::array<MAX_VOLUMES, intersection_t> intersected_volumes_dn_;

    int vol_dn;
    for (int i = 0; i < volumes_.size(); i++) {
      auto v = volumes_[i];
      auto hit_dn = v->intersects_with(dn);
      if (hit_dn.intersected) {
        intersected_volumes_dn_.emplace_back(i, hit_dn);
      }
    }

    util::heap_sort(intersected_volumes_dn_.begin(),
                    intersected_volumes_dn_.end(),
                    [&](const intersection_t& a, const intersection_t& b) {
                      return a.second.t_min < a.second.t_min;
                    });
    for (auto res : intersected_volumes_dn_) {
      F l_dn = res.second.t_max - res.second.t_min;
      F l_depth = model.deposition_depth(rng);
      if (l_depth <= l_dn) {
        hits++;

        response.detector2 = res.first;
        response.d2_entry = dn(res.second.t_min);
        response.d2_deposition = dn(res.second.t_min + l_depth);
        response.d2_exit = dn(res.second.t_max);

        break;
      }
    }

    return hits;
  }

 private:
  std::vector<PET3D::Full::Volume<F>*> volumes_;

  template <class RNG, class AcceptanceModel>
  bool find_first_interaction(const Ray& ray, HalfResponse response) {}
};
}
}

#endif  // SCANNER_H
