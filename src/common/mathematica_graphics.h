#pragma once

#include <iostream>
#include "2d/geometry/polygon.h"
#include "2d/barrel/detector_set.h"
#include "2d/geometry/line_segment.h"
#include "2d/geometry/pixel_grid.h"

namespace Common {

/// Outputs .m Mathematica file with given entities graphics
////
/// Creates .m Mathematica graphics for various classes.
template <typename FType> class MathematicaGraphics {
 public:
  using F = FType;

  MathematicaGraphics(std::ostream& out) : next_(false), out_(out) {
    out_ << "{\n";
  }
  ~MathematicaGraphics() { out_ << "}\n"; }

  template <std::size_t NumPoints>
  void add(const PET2D::Polygon<NumPoints, F>& polygon) {
    delimiter();
    out_ << "{Polygon[{\n";
    bool next = false;
    for (PET2D::Point<F> p : polygon) {
      delimiter(next) << "  " << pair(p.x, p.y);
    }
    out_ << "}]}";
  }

  template <class Detector, typename SType, std::size_t MaxDetectors>
  void add(const PET2D::Barrel::DetectorSet<Detector, SType, MaxDetectors>&
               scanner) {
    delimiter();
    out_ << "{";
    next_ = false;
    for (Detector detector : scanner) {
      add(detector);
    }
    out_ << "}";
    next_ = true;
  }

  template <class Detector, typename SType, std::size_t MaxDetectors>
  void add(
      const PET2D::Barrel::DetectorSet<Detector, SType, MaxDetectors>& scanner,
      const PET2D::Barrel::LOR<SType>& lor) {
    using F = typename Detector::F;
    delimiter();
    out_ << "{FaceForm[], EdgeForm[Black], MeshPrimitives[ConvexHullMesh[{\n";
    auto detector1 = scanner[lor.first];
    bool next1 = false;
    for (const PET2D::Point<F>& p : detector1) {
      delimiter(next1) << pair(p.x, p.y);
    }
    delimiter();
    auto detector2 = scanner[lor.second];
    bool next2 = false;
    for (const PET2D::Point<F>& p : detector2) {
      delimiter(next2) << pair(p.x, p.y);
    }
    out_ << "}], 2]}";
  }

  void add(const PET2D::LineSegment<F>& segment) {
    delimiter();
    out_ << "{Line[{";
    out_ << pair(segment.start.x, segment.start.y) << ", ";
    out_ << pair(segment.end.x, segment.end.y);
    out_ << "}]}";
  }

  void add_circle(const PET2D::Point<F>& center, F radius) {
    delimiter();
    out_ << "{Circle[";
    out_ << pair(center.x, center.y) << ", ";
    out_ << radius << "]}";
  }

  void add_circle(F radius) { add_circle(PET2D::Point<F>(0, 0), radius); }

  template <typename S>
  void add_pixel(const PET2D::PixelGrid<F, S>& grid,
                 const PET2D::Pixel<S> pixel) {
    auto ll = grid.lower_left_at(pixel);
    delimiter();
    out_ << "{FaceForm[], EdgeForm[Black], Polygon[{\n";
    out_ << pair(ll.x, ll.y) << ", ";
    out_ << pair(ll.x + grid.pixel_size, ll.y) << ", ";
    out_ << pair(ll.x + grid.pixel_size, ll.y + grid.pixel_size) << ", ";
    out_ << pair(ll.x, ll.y + grid.pixel_size) << "}]}";
  }

  void add(const PET2D::Point<F>& p) {
    delimiter();
    out_ << "{Red,Point[" << pair(p.x, p.y) << "]}";
  }

 private:
  std::string pair(F first, F second) {
    std::string result = "{";
    result += number(first) + ", " + number(second) + "}";
    return result;
  }

  std::string number(F number) {
    char number_char[64];
    sprintf(number_char, "%.12g", number);
    std::string number_str(number_char);
    auto i = number_str.find("e");
    if (i != std::string::npos) {
      number_str.erase(i, 1);
      number_str.insert(i, "*10^(");
      number_str.push_back(')');
    }
    return number_str;
  }

  std::ostream& delimiter(bool& next) {
    if (next) {
      out_ << ",\n";
    } else {
      next = true;
    }
    return out_;
  }

  std::ostream& delimiter() { return delimiter(next_); }

  bool next_;
  std::ostream& out_;
};

}  // Common
