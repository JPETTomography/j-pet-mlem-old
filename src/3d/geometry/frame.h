#pragma once

#include "3d/geometry/point.h"
#include "3d/geometry/vector.h"

template <typename FType> class Frame {
 public:
  using F = FType;
  using Point = PET3D::Point<F>;
  using Vector = PET3D::Vector<F>;
  using Point2D = PET2D::Point<F>;
  using Vector2D = PET2D::Vector<F>;

  Frame(const Point& p1, const Point& p2) : v_z_(0, 0, 1) {
    center_ = Point(0.5 * (p1.x + p2.x), 0.5 * (p1.y + p2.y), 0.0);
    Vector diff = p1 - p2;
    v_up_ = Vector(diff.x, diff.y, 0).normalized();
    v_normal_ = cross(v_up_, v_z_);
  }

  Frame(const Point2D& p1, const Point2D& p2) : v_z_(0, 0, 1) {
    center_ = Point(0.5 * (p1.x + p2.x), 0.5 * (p1.y + p2.y), 0.0);
    Vector2D diff = p1 - p2;
    v_up_ = Vector(diff.x, diff.y, 0).normalized();
    v_normal_ = cross(v_up_, v_z_);
  }

  Point2D project(const Point& p) {
    Vector v = p - center_;

    return Point2D(dot(v, v_z_), dot(v, v_up_));
  }

  F distance(const Point& p) {
    Vector v = p - center_;
    return dot(v, v_normal_);
  }

 private:
  Vector v_up_;
  Vector v_z_;
  Vector v_normal_;
  Point center_;
};
