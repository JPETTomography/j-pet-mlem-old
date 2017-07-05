#ifndef GATE_VOLUME_H
#define GATE_VOLUME_H

#include <list>

#include "2d/geometry/vector.h"

namespace Gate {
namespace D2 {

template <typename F> class Repeater {};

template <typename F> class Linear : public Repeater<F> {};
template <typename F> class Circular : public Repeater<F> {};

template <typename F> class Volume {
  using Vector = PET2D::Vector<F>;
  using VolumeList = std::list<Volume*>;

 public:
  Volume() : is_sd_(false) {}

  bool is_sd() const { return is_sd_; }

  void attach_daughter(Volume* daughter) { daughters_.push_back(daughter); }
  void attach_crystal_sd() { is_sd_ = true; }

  typename VolumeList::const_iterator daughters() const {
    return daughters_.begin();
  }
  typename VolumeList::const_iterator daughters_end() const {
    return daughters_.end();
  }

 private:
  // Daughters
  VolumeList daughters_;
  // Repeaters

  // Material
  // Translation
  // Rotation

  bool is_sd_;
};

template <typename F> class Box : public Volume<F> { public: };
}
}

#endif  // GATE_VOLUME_H
