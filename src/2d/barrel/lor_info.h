#ifndef LOR_INFO
#define LOR_INFO

#include "2d/barrel/lor.h"
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
  using PixelInfoContainer = std::vector<PixelInfo>;

  LorInfo(SType n_detectors)
      : n_detectors(n_detectors),
        max_index((int(n_detectors - 1) * (n_detectors)) / 2 + n_detectors - 2),
        lor_info_(max_index + 1) {}

  PixelInfoContainer& operator[](const LOR<SType>& lor) {
    return lor_info_[lor.index()];
  }

 public:
  SType n_detectors;
  int max_index;

 private:
  std::vector<PixelInfoContainer> lor_info_;
};
}
}
#endif  // LOR_INFO
