#ifndef LOR_INFO
#define LOR_INFO

#include "2d/barrel/lor.h"
#include "2d/geometry/pixel_grid.h"
#include "2d/geometry/point.h"
#include "2d/geometry/line_segment.h"

namespace PET2D {
namespace Barrel {

template <typename FType, typename SType> class LorPixelnfo {
 public:
  struct PixelInfo {
    PET2D::Pixel<SType> pixel;
    FType t;
    FType distance;
    FType fill;
  };

  using F = FType;
  using S = SType;
  using LOR = PET2D::Barrel::LOR<SType>;
  using PixelGrid = PET2D::PixelGrid<FType, SType>;
  using PixelInfoContainer = std::vector<PixelInfo>;
  using Point = PET2D::Point<F>;

  struct LorInfo {
    SType d1, d2;
    FType width;
    LineSegment<F> segment;
    PixelInfoContainer pixels;
  };

  LorPixelnfo(SType n_detectors, const PixelGrid& grid)
      : n_detectors(n_detectors),
        max_index((int(n_detectors - 1) * (n_detectors)) / 2 + n_detectors - 2),
        grid(grid),
        lor_info_(max_index + 1) {}

  LorInfo& operator[](const LOR& lor) { return lor_info_[lor.index()]; }

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
    auto segment = LineSegment<F>(Point(coords[2], coords[3]),
                                  Point(coords[0], coords[1]));
    F width;
    in.read((char*)&width, sizeof(F));
    int n_pixels;
    in.read((char*)&n_pixels, sizeof(S));
    if (in) {
      LOR lor(lor_desc[0], lor_desc[1]);
      lor_info_[lor.index()].width = width;
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

  const SType n_detectors;
  const int max_index;
  const PixelGrid grid;

 private:
  std::vector<LorInfo> lor_info_;
};
}
}
#endif  // LOR_INFO
