#ifndef GAUSS_SCANNER
#define GAUSS_SCANNER

namespace PET2D {
namespace Toy {

template <typename F> class GaussScanner {
 public:
  struct Response {
    F x, y;
  };

  using FullResponse = Response;

  GaussScanner(F s_z, F s_dl){};
};
}
}

#endif  // GAUSS_SCANNER
