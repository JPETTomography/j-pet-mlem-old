#ifndef GATE_VOLUME_H
#define GATE_VOLUME_H

#include <list>

#include "2d/geometry/vector.h"

namespace Gate {
namespace D2 {

template <typename F> class Repeater {};

template <typename F> class Linear : public Repeater<F> {};
template <typename F> class Circular : public Repeater<F> {};

template <typename FType> class Volume {
 public:
  using F = FType;
  using Vector = PET2D::Vector<F>;
  using VolumeList = std::list<Volume*>;

  Volume() : is_sd_(false) {}

  bool is_sd() const { return is_sd_; }

  typename VolumeList::const_iterator daughters() const {
    return daughters_.begin();
  }
  typename VolumeList::const_iterator daughters_end() const {
    return daughters_.end();
  }
  Vector translation() const { return translation_; }

  void attach_daughter(Volume* daughter) { daughters_.push_back(daughter); }
  void attach_crystal_sd() { is_sd_ = true; }

  void set_translation(Vector tr) { translation_ = tr; }

  virtual ~Volume() {}

 private:
  // Daughters
  VolumeList daughters_;
  // Repeaters

  // Material
  // Translation
  Vector translation_;
  // Rotation
  F angle_;

  bool is_sd_;
};

template <typename FType> class Box : public Volume<FType> {
 public:
  using F = FType;
  Box(F lX, F lY) : lengthX(lX), lengthY(lY) {}
  const F lengthX;
  const F lengthY;
};
}
}

#endif  // GATE_VOLUME_H
