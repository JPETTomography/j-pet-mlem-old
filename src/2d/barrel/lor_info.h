#pragma once

#include "2d/barrel/lor.h"
#include "2d/geometry/pixel_grid.h"
#include "2d/geometry/point.h"
#include "2d/geometry/line_segment.h"

namespace PET2D {
namespace Barrel {

template <typename FType, typename SType> struct LORInfo {
  using F = FType;
  using S = SType;
  using Pixel = PET2D::Pixel<S>;
  using LOR = PET2D::Barrel::LOR<S>;
  using LineSegment = PET2D::LineSegment<F>;

  struct PixelInfo {
    Pixel pixel;
    F t;
    F distance;
    F fill;
    F weight;
  };

  using PixelInfoList = std::vector<PixelInfo>;

  LOR lor;
  LineSegment segment;
  F width;
  PixelInfoList pixels;

  LORInfo() = default;

  LORInfo(const LOR& lor, const LineSegment& segment, const F width)
      : lor(lor), segment(segment), width(width) {}

  LORInfo(std::istream& in) : lor(in), segment(in), width(util::read<F>(in)) {
    auto n_pixels = util::read<int>(in);
    if (in && n_pixels) {
      pixels.resize(n_pixels);
      in.read((char*)&pixels[0], sizeof(PixelInfo) * n_pixels);
    }
  }

  void sort() {
    std::sort(pixels.begin(),
              pixels.end(),
              [](const PixelInfo& a, const PixelInfo& b) { return a.t < b.t; });
  }

  void push_back(const PixelInfo& pixel_info) { pixels.push_back(pixel_info); }
};

}  // Barrel
}  // PET2D
