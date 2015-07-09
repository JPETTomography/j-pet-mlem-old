#pragma once

#include "lor_info.h"

namespace PET2D {
namespace Barrel {

template <typename FType, typename SType>
class Geometry : public std::vector<LORInfo<FType, SType>> {
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

  Geometry(S n_detectors, const PixelGrid& grid)
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

  friend std::ostream& operator<<(std::ostream& out, const Geometry& lpi) {
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
