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

    HalfResponse response_up;
    int hits = find_first_interaction(rng, model, up, response_up);

    if (hits < 1)
      return hits;

    response.detector1 = response_up.detector;
    response.d1_entry = response_up.entry;
    response.d1_deposition = response_up.deposition;
    response.d1_exit = response_up.exit;

    HalfResponse response_dn;
    hits += find_first_interaction(rng, model, dn, response_dn);

    if (hits < 2)
      return hits;

    response.detector2 = response_dn.detector;
    response.d2_entry = response_dn.entry;
    response.d2_deposition = response_dn.deposition;
    response.d2_exit = response_dn.exit;

    return hits;
  }

 private:
  std::vector<PET3D::Full::Volume<F>*> volumes_;

  template <class RNG, class AcceptanceModel>
  int find_first_interaction(RNG& rng,
                             AcceptanceModel& model,
                             const Ray& ray,
                             HalfResponse& response) {
    using intersection_t = std::pair<int, ray_tracing::intersection_result<F>>;
    util::array<MAX_VOLUMES, intersection_t> intersected_volumes_;
    int hits = 0;

    for (int i = 0; i < volumes_.size(); i++) {
      auto v = volumes_[i];
      auto hit = v->intersects_with(ray);
      if (hit.intersected) {
        intersected_volumes_.emplace_back(i, hit);
      }
    }

    util::heap_sort(intersected_volumes_.begin(),
                    intersected_volumes_.end(),
                    [&](const intersection_t& a, const intersection_t& b) {
                      return a.second.t_min < b.second.t_min;
                    });
    for (auto res : intersected_volumes_) {
      F l = res.second.t_max - res.second.t_min;
      F l_depth = model.deposition_depth(rng);
      if (l_depth < l) {
        hits++;
        response.detector = res.first;
        response.entry = ray(res.second.t_min);
        response.deposition = ray(res.second.t_min + l_depth);
        response.exit = ray(res.second.t_max);
        break;
      }
    }
    return hits;
  }
};
}
}

#endif  // SCANNER_H
