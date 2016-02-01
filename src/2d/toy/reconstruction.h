#ifndef RECONSTRUCTION
#define RECONSTRUCTION

#include <vector>
namespace PET2D {
namespace Toy {

template <typename K, typename R> class Reconstruction {
  using Kernel = K;
  using Response = R;
  using F = Kernel::F;

  Event response(int i) const { return responses_[i]; }

 private:
  std::vector<Response> responses_;
};
}
}
#endif  // RECONSTRUCTION
