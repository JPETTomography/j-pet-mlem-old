#pragma once

#if !__CUDACC__
#include "util/bstream.h"
#endif

#include "2d/barrel/lor.h"
#include "2d/geometry/pixel_grid.h"
#include "2d/geometry/point.h"
#include "2d/geometry/line_segment.h"

namespace PET2D {
namespace Barrel {

/// Keeps extended information about LOR

/// Keeps extended information about LOR such as pixels belonging to the LOR
/// together with their description using PixelInfo.
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

#if !__CUDACC__
  // Construct LOR info from stream
  LORInfo(std::istream& in) : lor(in), segment(in), width(util::read<F>(in)) {
    size_t n_pixels;
    in >> n_pixels;
    if (in && n_pixels) {
      pixels.resize(n_pixels);
      in.read((char*)&pixels[0], sizeof(PixelInfo) * n_pixels);
    }
  }

  // Construct LOR info from binary stream
  LORInfo(util::ibstream& in) : lor(in), segment(in), width(in.read<F>()) {
    size_t n_pixels;
    in >> n_pixels;
    if (in && n_pixels) {
      pixels.resize(n_pixels);
      in.read((char*)&pixels[0], sizeof(PixelInfo) * n_pixels);
    }
  }

  friend util::obstream& operator<<(util::obstream& out,
                                    const LORInfo& lor_info) {
    out << lor_info.lor;
    out << lor_info.segment;
    out << lor_info.width;
    out << lor_info.pixels.size();
    out << lor_info.pixels;
    return out;
  }
#endif

  void sort() {
    std::sort(pixels.begin(),
              pixels.end(),
              [](const PixelInfo& a, const PixelInfo& b) { return a.t < b.t; });
  }

  void push_back(const PixelInfo& pixel_info) { pixels.push_back(pixel_info); }
};

}  // Barrel
}  // PET2D
