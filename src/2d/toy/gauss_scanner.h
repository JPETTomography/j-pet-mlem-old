#ifndef GAUSS_SCANNER
#define GAUSS_SCANNER

#include <random>

#include <2d/geometry/event.h>

namespace PET2D {
namespace Toy {

template <typename F> class GaussScanner {
 public:
  struct Response {
    F x, y;
  };

  using FullResponse = Response;
  using Event = PET2D::Event<F>;

  GaussScanner(F s_x, F s_y) : dist_x(0, s_x), dist_y(0,s_y){};
  template <class RNG>
  short exact_detect(RNG& rng,
                     const Event& event,
                     FullResponse& response) const {
    response.x = event.origin.x;
    response.y = event.origin.y;
    return 2;
  }

  template <class RNG, class AcceptanceModel>
  short exact_detect(RNG& rng,                ///< random number generator
                     AcceptanceModel& model,  ///< (unused)
                     const Event& event,      ///< event to be detected
                     Response& response       ///< scanner response (LOR+length)
                     ) const {
    (void)model;  // mark as unused
    return exact_detect(rng, event, response);
  }

  template <class RNG>
  Response response_w_error(RNG& rng, const FullResponse& full_response) {
    Response response = full_response;
    response.x += dist_x(rng);
    response.y += dist_y(rng);
    return response;
  }

  Response response_wo_error(const FullResponse& full_response) {
    return full_response;
  }

#if !__CUDACC__
  friend std::ostream& operator<<(std::ostream& out, const Response& response) {
    out << response.x << ' ' << response.y;
    return out;
  }
#endif

 private:
  std::normal_distribution<F> dist_x;
  std::normal_distribution<F> dist_y;
};
}
}

#endif  // GAUSS_SCANNER
