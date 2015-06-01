#ifndef LOR_INFO
#define LOR_INFO

#include "2d/barrel/lor.h"
#include "2d/geometry/pixel_grid.h"

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
  using LOR = PET2D::Barrel::LOR<SType>;
  using PixelGrid = PET2D::PixelGrid<FType, SType>;
  using PixelInfoContainer = std::vector<PixelInfo>;

  LorPixelnfo(SType n_detectors, const PixelGrid& grid)
      : n_detectors(n_detectors),
        max_index((int(n_detectors - 1) * (n_detectors)) / 2 + n_detectors - 2),
        grid(grid),
        lor_info_(max_index + 1) {}

  PixelInfoContainer& operator[](const LOR& lor) {
    return lor_info_[lor.index()];
  }

  std::istream& read(std::istream& in) {
    while (read_lor_info(in))
      ;
    return in;
  }

  std::istream& read_lor_info(std::istream& in) {
    int lor_desc[3];
    in.read((char*)lor_desc, 3 * sizeof(int));

    if (in) {
      LOR lor(lor_desc[0], lor_desc[1]);
      lor_info_[lor.index()].resize(lor_desc[2]);
      in.read((char*)&lor_info_[lor.index()][0],
              sizeof(PixelInfo) * lor_desc[2]);
    }
    return in;
  }

  void print(std::ostream& out) {
    for (int d1 = 0; d1 < n_detectors; ++d1) {
      for (int d2 = 0; d2 < d1; ++d2) {
        LOR lor(d1, d2);
        auto index = lor.index();
        if (lor_info_[index].size() > 0) {
          out << d1 << " " << d2 << "\n";
          for (PixelInfo& info : lor_info_[index]) {
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
  std::vector<PixelInfoContainer> lor_info_;
};
}
}
#endif  // LOR_INFO
