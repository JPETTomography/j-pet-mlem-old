#ifndef GATE_VOLUME_H
#define GATE_VOLUME_H

#include "2d/geometry/vector.h"

namespace Gate {
namespace D2 {

template <typename F> class Repeater {};

template <typename F> class Linear : public Repeater<F> {};
template <typename F> class Circular : public Repeater<F> {};

template <typename F> class Volume {
  using Vector = PET2D::Vector<F>;

 public:
  // Repeaters
  // Daughters
  // Material
  // Translation
  // Rotation
};

template <typename F> class Box : public Volume<F> { public: };
}
}

#endif  // GATE_VOLUME_H
