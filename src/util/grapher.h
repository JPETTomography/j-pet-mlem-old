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
  ~Graphics() { out_ << "}]\n"; }

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
      out_ << sep << "{" << p.x << "," << p.y << "}\n";
      sep = ",";
    }
    auto detector2 = scanner[lor.second];
    for (const PET2D::Point<F>& p : detector2) {
      out_ << sep << "{" << p.x << "," << p.y << "}\n";
      sep = ",";
    }
    out_ << "}],2]}\n";
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
