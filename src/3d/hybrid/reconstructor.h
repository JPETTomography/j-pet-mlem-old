#ifndef RECONSTRUCTOR
#define RECONSTRUCTOR

#include "2d/barrel/lor_info.h"

template <typename Scanner> class Reconstructor {
 public:
  using F = typename Scanner::F;
  using S = typename Scanner::S;
  using LorPixelInfo = PET2D::Barrel::LorPixelnfo<F, S>;
  using Response = typename Scanner::Response;
  using LOR = PET2D::Barrel::LOR<S>;

  struct Event {
    LOR lor;
    F up;
    F right;
    F tan;
  };

  Reconstructor(const Scanner& scanner, const LorPixelInfo& lor_pixel_info)
      : scanner_(scanner), lor_pixel_info_(lor_pixel_info) {}

 private:
  const Scanner& scanner_;
  const LorPixelInfo& lor_pixel_info_;
  std::vector<Response> responses_;
};

#endif  // RECONSTRUCTOR
