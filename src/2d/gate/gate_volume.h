#ifndef GATE_VOLUME_H
#define GATE_VOLUME_H

#include <list>

#include "2d/geometry/vector.h"

namespace Gate {
namespace D2 {

template <typename FType> class Repeater {
 public:
  using F = FType;
  using Vector = PET2D::Vector<F>;
};

template <typename FType> class Linear : public Repeater<FType> {
 public:
  using F = FType;
  using Vector = typename Repeater<F>::Vector;
  Linear(int n, const Vector& v){};
};
template <typename FType> class Circular : public Repeater<FType> {};

template <typename FType> class Volume {
 public:
  using F = FType;
  using Vector = PET2D::Vector<F>;
  using VolumeList = std::list<Volume*>;

  Volume() : is_sd_(false), translation_(0, 0), angle_(0) {}

  bool is_sd() const { return is_sd_; }

  typename VolumeList::const_iterator daughters() const {
    return daughters_.begin();
  }
  typename VolumeList::const_iterator daughters_end() const {
    return daughters_.end();
  }
  Vector translation() const { return translation_; }
  F rotation() const { return angle_; }

  void attach_daughter(Volume* daughter) { daughters_.push_back(daughter); }
  void attach_crystal_sd() { is_sd_ = true; }
  void attach_repeater(Repeater<F>* repeater){};
  void set_translation(Vector tr) { translation_ = tr; }
  void set_rotation(F a) { angle_ = a; }

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
