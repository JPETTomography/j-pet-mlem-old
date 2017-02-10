#ifndef SCANNER_H
#define SCANNER_H

#include <vector>

#include "3d/full/volume.h"
#include "3d/geometry/event.h"
#include "3d/geometry/ray.h"

namespace PET3D {
namespace Full {

template <typename F, typename S> class Scanner {
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

  void add_volume(Volume<F>* vol) {
    volumes_.push_back(vol);
    intersected_volumes_.reserve(volumes_.size());
  };

  template <class RNG, class AcceptanceModel>
  S exact_detect(RNG& rng,                ///< random number generator
                 AcceptanceModel& model,  ///< acceptance model
                 const Event& event,      ///< event to be detected
                 FullResponse& response   ///< response (LOR+zu+zd+dl)
                 ) {
    return 0;
    ray_tracing::Ray<F> up(event.origin, event.direction);
    ray_tracing::Ray<F> dn(event.origin, -event.direction);
    S intersected_up = 0;

    for (int i = 0; i < volumes_.size(); i++) {
      auto v = volumes_[i];
      auto hit_up = v->intersects_with(up);
      if (hit_up.intersected) {
        intersected_up++;
      }
    }
    //    F l_up = hit_up.t_max - hit_up.t_min;
    //    F l_depth_up = model.deposition_depth(rng);
    //    if (l_depth_up > l_up)
    //      return 0;
  }

 private:
  std::vector<PET3D::Full::Volume<F>*> volumes_;
  std::vector<std::pair<int, ray_tracing::intersection_result<F>>>
      intersected_volumes_;
};
}
}

#endif  // SCANNER_H
