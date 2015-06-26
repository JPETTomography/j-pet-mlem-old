#ifndef LM_RECONSTRUCTION
#define LM_RECONSTRUCTION

#include "2d/barrel/lors_pixels_info.h"
#include "2d/barrel/detector_set.h"
#include "2d/geometry/point.h"
#include "2d/barrel/lor.h"

namespace PET2D {
namespace Barrel {

template <typename FType, typename SType> class LMReconstruction {
  using Detector = PET2D::Barrel::SquareDetector<FType>;
  using Scanner2D = PET2D::Barrel::DetectorSet<Detector, 192, SType>;

 public:
  using F = FType;
  using S = SType;
  using Response = typename Scanner2D::Response;
  using Point = PET2D::Point<F>;
  using LOR = PET2D::Barrel::LOR<S>;
  using LORsPixelsInfo = PET2D::Barrel::LORsPixelsInfo<F, S>;
  using PixelInfo = typename LORsPixelsInfo::PixelInfo;
  using PixelConstIterator =
      typename LORsPixelsInfo::PixelInfoContainer::const_iterator;

  struct BarrelEvent {
    LOR lor;
    Point p;
    PixelConstIterator first_pixel;
    PixelConstIterator last_pixel;

    F gauss_norm;
    F inv_sigma2;
  };

  LMReconstruction(const LORsPixelsInfo& lor_pixel_info, F sigma)
      : lor_pixel_info(lor_pixel_info) {}

  int fscanf_responses(std::istream& in) {
    int count = 0;
    while (true) {
      auto response = fscanf_response(in);
      if (!in)
        break;
      count++;
      BarrelEvent event;
      event.lor = response.lor;

      auto segment = lor_pixel_info[response.lor].segment;
      F t = 0.5 - response.dl / (2 * segment->length);
      event.p = PET2D::interpolate(t, segment->start, segment->end);

      PixelInfo pix_info_up, pix_info_dn;
      pix_info_up.t = t + 3 * sigma_;
      pix_info_dn.t = t - 3 *sigma_;
      event.last_pixel =
          std::upper_bound(lor_pixel_info[event.lor].pixels.begin(),
                           lor_pixel_info[event.lor].pixels.end(),
                           pix_info_up,
                           [](const PixelInfo& a, const PixelInfo& b)
                               -> bool { return a.t < b.t; });

      event.first_pixel =
          std::lower_bound(lor_pixel_info[event.lor].pixels.begin(),
                           lor_pixel_info[event.lor].pixels.end(),
                           pix_info_dn,
                           [](const PixelInfo& a, const PixelInfo& b)
                               -> bool { return a.t < b.t; });
    }
    return count;
  }

  const LORsPixelsInfo& lor_pixel_info;

 private:
  Response fscanf_response(std::istream& in) {
    S d1, d2;
    Response response;
    in >> d1 >> d2 >> response.tof_position >> response.dl;
    response.lor = LOR(d1, d2);
    return response;
  }

  std::vector<BarrelEvent> events_;
  F sigma_;
};
}
}
#endif  // LM_RECONSTRUCTION
