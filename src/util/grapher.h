#pragma once

#include <iostream>
#include "2d/geometry/polygon.h"
#include "2d/barrel/detector_set.h"

template <typename FType> class Graphics {
 public:
  using F = FType;

  Graphics(std::ostream& out) : first_(true), out_(out) {
    out_ << "Graphics[{\n";
  };
  ~Graphics() { out_ << "}]\n";}

  template <std::size_t NumPoints>
  void add(const PET2D::Polygon<NumPoints, F>& polygon) {
    add();
    out_ << "{Polygon[{";
    std::string sep("");
    for (PET2D::Point<F> p : polygon) {
      out_ << sep << "{" << p.x << "," << p.y << "}\n";
      sep = ",";
    }
    out_ << "}]}\n";
  }


 private:
  void add() {
    if (first_) {
      first_ = false;
    } else
      out_ << ',';
  }

  bool first_;
  std::ostream& out_;
};
