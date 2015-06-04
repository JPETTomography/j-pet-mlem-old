#pragma once

#include <iostream>
#include "2d/geometry/polygon.h"
#include "2d/barrel/detector_set.h"
#include "2d/geometry/line_segment.h"

template <typename FType> class Graphics {
 public:
  using F = FType;

  Graphics(std::ostream& out) : first_(true), out_(out) {
    out_ << "Graphics[{\n";
  };
  ~Graphics() { out_ << "}]\n"; }

  template <std::size_t NumPoints>
  void add(const PET2D::Polygon<NumPoints, F>& polygon) {
    add();
    out_ << "{Polygon[{";
    std::string sep("");
    for (PET2D::Point<F> p : polygon) {
      out_ << sep << pair(p.x, p.y) << "\n";
      sep = ",";
    }
    out_ << "}]}\n";
  }

  template <typename Detector, std::size_t MaxDet, typename SType>
  void add(const PET2D::Barrel::DetectorSet<Detector, MaxDet, SType>& scanner) {
    add();
    out_ << "{\n";
    first_ = true;
    for (Detector detector : scanner) {
      add(detector);
    }
    out_ << "}\n";
  }

  template <typename Detector, std::size_t MaxDet, typename SType>
  void add(const PET2D::Barrel::DetectorSet<Detector, MaxDet, SType>& scanner,
           const PET2D::Barrel::LOR<SType>& lor) {
    using F = typename Detector::F;
    add();
    out_ << "{FaceForm[], EdgeForm[Black], MeshPrimitives[ConvexHullMesh[{\n";
    std::string sep = "";
    auto detector1 = scanner[lor.first];
    for (const PET2D::Point<F>& p : detector1) {
      out_ << sep << pair(p.x, p.y) << "\n";
      sep = ",";
    }
    auto detector2 = scanner[lor.second];
    for (const PET2D::Point<F>& p : detector2) {
      out_ << sep << pair(p.x, p.y) << "\n";
      sep = ",";
    }
    out_ << "}],2]}\n";
  }

  void add(const PET2D::LineSegment<F>& segment) {
    add();
    out_ << "{Line[{";
    out_ << pair(segment.start.x, segment.start.y) << ",";
    out_ << pair(segment.end.x, segment.end.y);
    out_ << "}]}";
  }

  void addCircle(const PET2D::Point<F>& center, F radius) {
    add();
    out_ << "{Circle[";
    out_ << pair(center.x, center.y) << ",";
    out_ << radius << "]}\n";
  }

  void addCircle(F radius) { addCircle(PET2D::Point<F>(0, 0), radius); }

 private:
  std::string pair(F first, F second) {
    std::string result = "{";
    result += number(first) + "," + number(second) + "}";
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

  void add() {
    if (first_) {
      first_ = false;
    } else
      out_ << ',';
  }

  bool first_;
  std::ostream& out_;
};
