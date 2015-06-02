#ifndef RECONSTRUCTOR
#define RECONSTRUCTOR

#include <vector>

#include "2d/barrel/lor_info.h"
#include "2d/strip/event.h"

template <typename Scanner> class Reconstructor {
 public:
  using F = typename Scanner::F;
  using S = typename Scanner::S;
  using LorPixelInfo = PET2D::Barrel::LorPixelnfo<F, S>;
  using Response = typename Scanner::Response;
  using LOR = PET2D::Barrel::LOR<S>;
  using StripEvent = PET2D::Strip::Event<F>;

  struct FrameEvent {
    LOR lor;
    F up;
    F right;
    F tan;
  };

  Reconstructor(const Scanner& scanner, const LorPixelInfo& lor_pixel_info)
      : scanner_(scanner), lor_pixel_info_(lor_pixel_info) {}

  FrameEvent translate_to_frame(const Response& response) {
    FrameEvent event;
    event.lor = response.lor;

    auto R = lor_pixel_info_[event.lor].segment.length;
    StripEvent strip_event(response.z_up, response.z_dn, response.dl);
    F tan, y, z;
    strip_event.transform(R, event.tan, event.up, event.right);
    return event;
  }

  Response fscanf_responses(std::istream& in) {
    auto response = fscanf_response(in);
    while (in) {
      events_.push_back(translate_to_frame(response));
      response = fscanf_response(in);
    }
  }

  int n_events() const { return events_.size(); }
  FrameEvent frame_event(int i) const { return events_[i]; }

 private:
  Response fscanf_response(std::istream& in) {
    S d1, d2;
    F z_up, z_dn, dl;
    Response response;
    in >> d1 >> d2 >> response.z_up >> response.z_dn >> response.dl;
    response.lor = LOR(d1, d2);
    return response;
  }

  const Scanner& scanner_;
  const LorPixelInfo& lor_pixel_info_;
  std::vector<FrameEvent> events_;
};

#endif  // RECONSTRUCTOR
