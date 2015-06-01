#ifndef LOR_INFO
#define LOR_INFO

#include "2d/barrel/lor.h"
#include "2d/geometry/pixel_grid.h"

namespace PET2D {
namespace Barrel {

template <typename FType, typename SType> class LorInfo {
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

  LorInfo(SType n_detectors, const PixelGrid& grid)
      : n_detectors(n_detectors),
        max_index((int(n_detectors - 1) * (n_detectors)) / 2 + n_detectors - 2),
        lor_info_(max_index + 1),
        grid(grid) {}

  PixelInfoContainer& operator[](const LOR& lor) {
    return lor_info_[lor.index()];
  }

  std::ifstream& read(std::ifstream& in) {
    while (read_lor_info(in))
      ;
    return in;
  }

  std::ifstream& read_lor_info(std::ifstream& in) {
    int lor_desc[3];
    in.read((char*)lor_desc, 3 * sizeof(int));
    if (in) {
      LOR lor(lor_desc[0], lor_desc[1]);
      lor_info_[lor].reserve(lor_desc[2]);
      in.read((char*)&lor_info_[lor][0], sizeof(PixelInfo) * lor_desc[2]);
    }
    return in;
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
