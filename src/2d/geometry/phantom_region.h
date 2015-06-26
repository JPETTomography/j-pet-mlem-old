#ifndef PHANTOM_REGION
#define PHANTOM_REGION

namespace PET2D {

// template <typename F> int sgn(F val) { return (0 < val) - (val < 0); }

/// Virtual phantom region made of ellipse and intensity
template <typename FType, typename RNG> struct PhantomRegion {
  using F = FType;
  using Point = PET2D::Point<F>;

  PhantomRegion(F intensity) : intensity(intensity) {}

  virtual bool contains(Point p) const = 0;
  virtual F weight() const = 0;
  virtual Point random_point(RNG&) = 0;

  const F intensity;
};

template <typename Shape, typename RNG>
class ShapePhantomRegion : public PhantomRegion<typename Shape::F, RNG> {
 public:
  using F = typename Shape::F;

  ShapePhantomRegion(const Shape& shape, F intensity)
      : PhantomRegion<F, RNG>(intensity),
        shape(shape),
        weight_(intensity * shape.area) {}

  bool contains(Point<F> p) const { return shape.contains(p); }

  const Shape shape;

  F weight() const { return weight_; }

 private:
  const F weight_;
};

template <typename FType, typename RNG>
class EllipticalPhantomRegion
    : public ShapePhantomRegion<PET2D::Ellipse<FType>, RNG> {
 public:
  using F = FType;
  using Ellipse = PET2D::Ellipse<F>;
  using Point = PET2D::Point<F>;

  EllipticalPhantomRegion(const Ellipse& ellipse, F intensity)
      : ShapePhantomRegion<Ellipse, RNG>(ellipse, intensity), gen_(ellipse) {}

  Point random_point(RNG& rng) { return gen_(rng); }

 private:
  EllipsePointGenerator<F> gen_;
};

template <typename FType, typename RNG>
class RectangularPhantomRegion
    : public ShapePhantomRegion<PET2D::Rectangle<FType>, RNG> {
 public:
  using F = FType;
  using Rectangle = PET2D::Rectangle<F>;
  using Point = PET2D::Point<F>;

  RectangularPhantomRegion(const Rectangle& rectangle, F intensity)
      : ShapePhantomRegion<Rectangle, RNG>(rectangle, intensity),
        gen_(rectangle) {}

  Point random_point(RNG& rng) { return gen_(rng); }

 private:
  PET2D::RectanglePointGenerator<F> gen_;
};
}


#endif  // PHANTOM_REGION
