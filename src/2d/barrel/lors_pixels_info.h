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

  struct PixelInfo {
    Pixel pixel;
    F t;
    F distance;
    F fill;
    F weight;
  };

  using PixelInfoList = std::vector<PixelInfo>;

  LOR lor;
  LineSegment<F> segment;
  F width;
  PixelInfoList pixels;

  LORInfo() = default;

  LORInfo(std::istream& in) : lor(in), segment(in), width(util::read<F>(in)) {
    auto n_pixels = util::read<int>(in);
    if (in && n_pixels) {
      pixels.resize(n_pixels);
      in.read((char*)&pixels[0], sizeof(PixelInfo) * n_pixels);
    }
  }
};

template <typename FType, typename SType> class LORsPixelsInfo {
 public:
  using F = FType;
  using S = SType;

  using Pixel = PET2D::Pixel<S>;
  using LOR = PET2D::Barrel::LOR<S>;
  using PixelGrid = PET2D::PixelGrid<F, S>;
  using Point = PET2D::Point<F>;
  using LORInfo = PET2D::Barrel::LORInfo<F, S>;
  using PixelInfo = typename LORInfo::PixelInfo;

  LORsPixelsInfo(S n_detectors, const PixelGrid& grid)
      : n_detectors(n_detectors),
        max_index((int(n_detectors - 1) * (n_detectors)) / 2 + n_detectors - 2),
        grid(grid),
        lor_info_(max_index + 1) {}
  LORInfo& operator[](const LOR& lor) { return lor_info_[lor.index()]; }

  const LORInfo& operator[](const LOR& lor) const {
    return lor_info_[lor.index()];
  }

  void push_back_pixel_info(const LOR& lor, const PixelInfo& pinfo) {
    lor_info_[lor.index()].pixels.push_back(pinfo);
  }

  void push_back_pixel(const LOR& lor, const Pixel& pixel, F weight) {

    auto& lor_info = lor_info_[lor.index()];
    auto center = grid.center_at(pixel.x, pixel.y);
    auto t = lor_info.segment.projection_scaled(center);
    PixelInfo pinfo;
    pinfo.pixel = pixel;
    pinfo.t = t;
    pinfo.weight = weight;

    push_back_pixel_info(lor, pinfo);
  }

  void sort(LORInfo& lor_info) {
    std::sort(lor_info.pixels.begin(),
              lor_info.pixels.end(),
              [](const PixelInfo& a, const PixelInfo& b) { return a.t < b.t; });
  }

  void sort() {
    for (auto& lor_info : lor_info_) {
      sort(lor_info);
    }
  }

  // assumes sorted pixels
  void merge_duplicates(LORInfo& lor_info) {
    auto& pixels = lor_info.pixels;
    for (auto& pinfo : pixels) {
    }
  }

  // Reading (binary)
  void read(std::istream& in) {
    while (in) {
      LORInfo lor_info(in);
      lor_info_[lor_info.lor.index()] = std::move(lor_info);
    }
  }

  friend std::ostream& operator<<(std::ostream& out,
                                  const LORsPixelsInfo& lpi) {
    for (int d1 = 0; d1 < lpi.n_detectors; ++d1) {
      for (int d2 = 0; d2 < d1; ++d2) {
        LOR lor(d1, d2);
        auto index = lor.index();
        if (lpi.lor_info_[index].pixels.size() > 0) {
          out << d1 << " " << d2 << " " << lpi.lor_info_[index].width
              << std::endl;
          for (PixelInfo& info : lpi.lor_info_[index].pixels) {
            out << info.pixel.x << " " << info.pixel.y << " " << info.t << " "
                << info.distance << " " << info.fill << std::endl;
          }
        }
      }
    }
    return out;
  }

  // Dirty hack
  void erase_pixel_info() {
    for (auto& lor_info : lor_info_) {
      lor_info.pixels.clear();
    }
  }

  const S n_detectors;
  const int max_index;
  const PixelGrid grid;

  typename std::vector<LORInfo>::const_iterator begin() const {
    return lor_info_.begin();
  }
  typename std::vector<LORInfo>::iterator begin() { return lor_info_.begin(); }
  typename std::vector<LORInfo>::const_iterator end() const {
    return lor_info_.end();
  }
  typename std::vector<LORInfo>::iterator end() { return lor_info_.end(); }

 private:
  std::vector<LORInfo> lor_info_;
};
}
}
