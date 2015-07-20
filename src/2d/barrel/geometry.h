#pragma once

#if !__CUDACC__
#include "util/bstream.h"
#include "util/read.h"
#endif

#include "lor_info.h"

namespace PET2D {
namespace Barrel {

#if !__CUDACC__
namespace {
template <typename FType> constexpr uint32_t magic() { return 0; }
template <> constexpr uint32_t magic<float>() { return "PETg"_4cc; }
template <> constexpr uint32_t magic<double>() { return "PETG"_4cc; }
}
#endif

/// Keeps extended information about barrel, including grid and all LOR info
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

#if !__CUDACC__
  /// Construct geometry from stream
  Geometry(std::istream& in) : Geometry(util::read<uint32_t>(in), in) {}

  /// Construct geometry from binary stream
  Geometry(util::ibstream& in) : Geometry(in.read<uint32_t>(), in) {}

  /// Serialize geometry into binary output stream
  friend util::obstream& operator<<(util::obstream& out,
                                    const Geometry& geometry) {
    out << magic<F>();
    out << geometry.n_detectors;
    out << geometry.grid;
    for (const auto& lor_info : geometry) {
      out << lor_info;
    }
    return out;
  }

 private:
  Geometry(S n_detectors, std::istream& in)
      : Base(((int(n_detectors) + 1) * (n_detectors)) / 2),
        n_detectors(n_detectors),
        grid(in) {
    while (in) {
      LORInfo lor_info(in);
      (*this)[lor_info.lor] = std::move(lor_info);
    }
  }
  Geometry(uint32_t file_magic, S n_detectors, util::ibstream& in)
      : Base(((int(n_detectors) + 1) * (n_detectors)) / 2),
        n_detectors(n_detectors),
        grid(in) {
    if (file_magic != magic<F>())
      throw("unknown binary geometry file");
    while (in) {
      LORInfo lor_info(in);
      (*this)[lor_info.lor] = std::move(lor_info);
    }
  }
  Geometry(uint32_t file_magic, util::ibstream& in)
      : Geometry(file_magic, in.read<S>(), in) {}

 public:
#endif

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
