#pragma once

#include "2d/barrel/lor.h"
#include "2d/geometry/pixel_grid.h"
#include "2d/geometry/point.h"
#include "2d/geometry/line_segment.h"

namespace PET2D {
namespace Barrel {

template <typename FType, typename SType> class LORsPixelsInfo {
 public:
  using F = FType;
  using S = SType;

  using Pixel = PET2D::Pixel<S>;
  using LOR = PET2D::Barrel::LOR<S>;
  using PixelGrid = PET2D::PixelGrid<F, S>;
  using Point = PET2D::Point<F>;

  struct PixelInfo {
    Pixel pixel;
    F t;
    F distance;
    F fill;
  };

  using PixelInfoContainer = std::vector<PixelInfo>;

  struct LORInfo {
    S d1, d2;
    F width;
    PixelInfoContainer pixels;
    LineSegment<F>* segment;
  };

  LORsPixelsInfo(S n_detectors, const PixelGrid& grid)
      : n_detectors(n_detectors),
        max_index((int(n_detectors - 1) * (n_detectors)) / 2 + n_detectors - 2),
        grid(grid),
        lor_info_(max_index + 1) {}

  LORInfo& operator[](const LOR& lor) { return lor_info_[lor.index()]; }

  const LORInfo& operator[](const LOR& lor) const {
    return lor_info_[lor.index()];
  }
  std::istream& read(std::istream& in) {
    while (read_lor_info(in))
      ;
    return in;
  }

  std::istream& read_lor_info(std::istream& in) {
    int lor_desc[2];
    in.read((char*)lor_desc, 2 * sizeof(S));
    F coords[4];
    in.read((char*)coords, 4 * sizeof(F));
    //    std::cout << coords[0] << " " << coords[1] << " " << coords[2] << " "
    //              << coords[3] << "\n";

    F width;
    in.read((char*)&width, sizeof(F));
    int n_pixels;
    in.read((char*)&n_pixels, sizeof(S));
    if (in) {
      LOR lor(lor_desc[0], lor_desc[1]);
      lor_info_[lor.index()].width = width;
      lor_info_[lor.index()].segment = new LineSegment<F>(
          Point(coords[2], coords[3]), Point(coords[0], coords[1]));
      lor_info_[lor.index()].pixels.resize(n_pixels);

      in.read((char*)&lor_info_[lor.index()].pixels[0],
              sizeof(PixelInfo) * n_pixels);
    }
    return in;
  }

  void print(std::ostream& out) {
    for (int d1 = 0; d1 < n_detectors; ++d1) {
      for (int d2 = 0; d2 < d1; ++d2) {
        LOR lor(d1, d2);
        auto index = lor.index();
        if (lor_info_[index].pixels.size() > 0) {
          out << d1 << " " << d2 << " " << lor_info_[index].width << "\n";
          for (PixelInfo& info : lor_info_[index].pixels) {
            out << info.pixel.x << " " << info.pixel.y << " " << info.t << " "
                << info.distance << " " << info.fill << "\n";
          }
        }
      }
    }
  }

  const S n_detectors;
  const int max_index;
  const PixelGrid grid;

 private:
  std::vector<LORInfo> lor_info_;
};
}
}
