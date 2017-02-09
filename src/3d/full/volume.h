#ifndef VOLUME_H
#define VOLUME_H

#include "3d/geometry/ray.h"

namespace PET3D {
namespace Full {

template <typename F> class Material {};

template <typename F> class Volume {
 public:
  virtual ray_tracing::intersection_result<F> intersects_with(
      const ray_tracing::Ray<F>& ray) = 0;
  bool is_detector;

 private:
  Material<F> material;
};

template <typename F> class BoxVolume : public Volume<F> {
  virtual ray_tracing::intersection_result<F> intersects_with(
      const ray_tracing::Ray<F>& ray) {
    return ray_tracing::intersect(ray, box);
  }

 private:
  ray_tracing::Box<F> box;
};
}
}
#endif  // VOLUME_H
