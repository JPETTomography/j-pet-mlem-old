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

    ray_tracing::Ray<F> up(event.origin, event.direction);
    ray_tracing::Ray<F> dn(event.origin, -event.direction);
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
      F l_depth_up = model.deposition_depth(rng);
      if (l_depth_up < l_up) {
        hits++;
        vol_up = res.first;
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
      F l_depth_dn = model.deposition_depth(rng);
      if (l_depth_dn <= l_dn) {
        hits++;
        vol_dn = res.first;
        break;
      }
    }

    return hits;
  }

 private:
  std::vector<PET3D::Full::Volume<F>*> volumes_;
};
}
}

#endif  // SCANNER_H
