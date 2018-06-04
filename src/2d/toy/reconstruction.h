#pragma once

#include <vector>

namespace PET2D {
namespace Toy {

template <typename K, typename R> class Reconstruction {

 public:
  using Kernel = K;
  using Response = R;
  using F = typename Kernel::F;

  Response response(int i) const { return responses_[i]; }

  template <typename StreamType> Reconstruction& operator<<(StreamType& in) {
    int i = 0;
    for (;;) {
      F x, y;
      in >> x >> y;
      if (!in)
        break;
      responses_.emplace_back(x, y);
      i++;
    }
    return *this;
  }

  size_t n_resp() const { return responses_.size(); }

 private:
  std::vector<Response> responses_;
};
}
}
