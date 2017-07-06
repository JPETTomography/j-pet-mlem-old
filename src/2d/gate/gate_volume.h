#ifndef GATE_VOLUME_H
#define GATE_VOLUME_H

#include <list>

#include "2d/geometry/vector.h"
#include "2d/geometry/transformation.h"

namespace Gate {
namespace D2 {

template <typename FType> class Volume;

template <typename FType> class Repeater {
 public:
  using F = FType;
  using Vector = PET2D::Vector<F>;
  using Transformation = PET2D::Transformation<F>;
  Repeater(int n) : n(n) {}
  const int n;
  virtual Transformation operator[](int i) = 0;
};

template <typename FType> class Linear : public Repeater<FType> {
 public:
  using F = FType;
  using Vector = typename Repeater<F>::Vector;
  using Transformation = PET2D::Transformation<F>;

  Linear(int n, const Vector& v)
      : Repeater<FType>(n), v(v), transformations_() {
    Vector center = v * (n - 1) * 0.5;
    for (int i = 0; i < n; i++) {
      transformations_.push_back(Transformation(0, v * i - center));
    }
  };
  Transformation operator[](int i) { return transformations_[i]; }
  const Vector v;

 public:
  std::vector<Transformation> transformations_;
};

template <typename FType> class Ring : public Repeater<FType> {
 public:
  using F = FType;
  using Vector = typename Repeater<F>::Vector;
  using Transformation = PET2D::Transformation<F>;

  Ring(int n, F angle, Vector center) : Repeater<FType>(n), center(center) {
    for (int i = 0; i < n; i++) {
      F d_angle = 2 * M_PI / n;
      transformations_.push_back(Transformation(center) *
                                 Transformation(i * d_angle) *
                                 Transformation(-center));
    }
  }

  Transformation operator[](int i) { return transformations_[i]; }

  const Vector center;

 public:
  std::vector<Transformation> transformations_;
};

template <typename FType> class Volume {
 public:
  using F = FType;
  using Vector = PET2D::Vector<F>;
  using VolumeList = std::list<Volume*>;
  using Transformation = PET2D::Transformation<F>;

  Volume() : is_sd_(false), transformation_(new Transformation()) {}

  bool is_sd() const { return is_sd_; }

  typename VolumeList::const_iterator daughters() const {
    return daughters_.begin();
  }
  typename VolumeList::const_iterator daughters_end() const {
    return daughters_.end();
  }

  Transformation transformation() const { return *transformation_; }
  Repeater<F>* repeater() const { return repeater_.get(); }

  void attach_daughter(Volume* daughter) { daughters_.push_back(daughter); }
  void attach_crystal_sd() { is_sd_ = true; }
  void attach_repeater(Repeater<F>* repeater) {
    repeater_ = std::unique_ptr<Repeater<F>>(repeater);
  };
  void detach_repeater() { repeater_.release(); }
  void set_translation(Vector tr) {
    transformation_ = std::unique_ptr<Transformation>(
        new Transformation(transformation_->rotation, tr));
  }
  void set_rotation(F a) {
    transformation_ = std::unique_ptr<Transformation>(
        new Transformation(a, transformation_->translation));
  }

  void set_placement(F a, Vector t) {
    transformation_ = std::unique_ptr<Transformation>(new Transformation(a, t));
  }

  void set_transformation(Transformation* t) {
    transformation_ = std::unique_ptr<Transformation>(t);
  }

  virtual ~Volume() {}

 private:
  VolumeList daughters_;
  std::unique_ptr<Repeater<F>> repeater_;
  // Material
  std::unique_ptr<Transformation> transformation_;

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
