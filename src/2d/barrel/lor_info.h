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

template <typename FType, typename SType>
class LORInfoList : public std::vector<LORInfo<FType, SType>> {
 public:
  using F = FType;
  using S = SType;

  using LORInfo = PET2D::Barrel::LORInfo<F, S>;
  using Base = std::vector<LORInfo>;
  using Pixel = PET2D::Pixel<S>;
  using LOR = PET2D::Barrel::LOR<S>;
  using PixelGrid = PET2D::PixelGrid<F, S>;
  using Point = PET2D::Point<F>;
  using PixelInfo = typename LORInfo::PixelInfo;

  LORInfoList(S n_detectors, const PixelGrid& grid)
      : Base(((int(n_detectors) + 1) * (n_detectors)) / 2),
        n_detectors(n_detectors),
        grid(grid) {}

  LORInfo& operator[](const LOR& lor) { return this->at(lor.index()); }

  const LORInfo& operator[](const LOR& lor) const {
    return this->at(lor.index());
  }

  void push_back_pixel_info(const LOR& lor, const PixelInfo& pixel_info) {
    (*this)[lor].pixels.push_back(pixel_info);
  }

  void push_back_pixel(const LOR& lor, const Pixel& pixel, F weight) {
    auto& lor_info = (*this)[lor];
    auto center = grid.center_at(pixel.x, pixel.y);
    auto t = lor_info.segment.projection_scaled(center);
    PixelInfo pixel_info;
    pixel_info.pixel = pixel;
    pixel_info.t = t;
    pixel_info.weight = weight;
    push_back_pixel_info(lor, pixel_info);
  }

  void sort_all() {
    for (auto& lor_info : *this) {
      lor_info.sort();
    }
  }

  void erase_pixel_info() {
    for (auto& lor_info : *this) {
      lor_info.pixels.resize(0);
    }
  }

  // Reading (binary)
  void read(std::istream& in) {
    while (in) {
      LORInfo lor_info(in);
      (*this)[lor_info.lor] = std::move(lor_info);
    }
  }

  friend std::ostream& operator<<(std::ostream& out, const LORInfoList& lpi) {
    for (int d1 = 0; d1 < lpi.n_detectors; ++d1) {
      for (int d2 = 0; d2 < d1; ++d2) {
        LOR lor(d1, d2);
        auto index = lor.index();
        if (lpi[index].pixels.size() > 0) {
          out << d1 << " " << d2 << " " << lpi[index].width << std::endl;
          for (PixelInfo& info : lpi[index].pixels) {
            out << info.pixel.x << " " << info.pixel.y << " " << info.t << " "
                << info.distance << " " << info.fill << std::endl;
          }
        }
      }
    }
    return out;
  }

  const S n_detectors;
  const PixelGrid grid;
};

}  // Barrel
}  // PET2D
