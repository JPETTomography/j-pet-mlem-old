#ifndef SCANNER_H
#define SCANNER_H

#include "3d/full/volume.h"
#include "3d/geometry/event.h"

namespace PET3D {
namespace Full {

template <typename F, typename S> class Scanner {
  using Event = PET3D::Event<F>;
  using Point = PET3D::Point<F>;

  struct FullResponse {
    S detector1, detector2;
    Point d1_entry, d1_exit, d1_deposition;
    Point d2_entry, d2_exit, d2_deposition;
    Point origin;

    FullResponse() = default;
  };

  void add_volume(const Volume<F>& vol);

  template <class RNG, class AcceptanceModel>
  short exact_detect(RNG& rng,                ///< random number generator
                     AcceptanceModel& model,  ///< acceptance model
                     const Event& event,      ///< event to be detected
                     FullResponse& response   ///< response (LOR+zu+zd+dl)
                     );
};
}
}

#endif  // SCANNER_H
