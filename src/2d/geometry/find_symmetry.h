#ifndef FIND_SYMMETRY_H
#define FIND_SYMMETRY_H
#include "2d/geometry/transformation.h"
#include "2d/barrel/symmetry_descriptor.h"

template <typename F, typename S>
PET2D::Transformation<F> symmetry_transformation(S symmetry) {
  using Transformation = PET2D::Transformation<F>;

  using Vector = typename Transformation::Vector;
  switch (symmetry) {
    case 0:
      return Transformation();
    case 1:
      return Transformation(M_PI, Vector(), true);
    case 2:
      return Transformation(0, Vector(), true);
    case 3:
      return Transformation(M_PI, Vector(), false);
    case 4:
      return Transformation(-M_PI / 2, Vector(), true);
    case 5:
      return Transformation(M_PI / 2, Vector(), false);
    case 6:
      return Transformation(-M_PI / 2, Vector(), false);
    case 7:
      return Transformation(M_PI / 2, Vector(), true);
  }
}

template <typename Scanner>
typename Scanner::S find_symmetric(Scanner scanner,
                                   typename Scanner::S symmetry,
                                   typename Scanner::S detector,
                                   typename Scanner::F epsilon) {
  using F = typename Scanner::F;
  using S = typename Scanner::S;

  auto t = symmetry_transformation<F>(symmetry);
  auto t_detector = scanner[detector].transformed(t);

  for (S d = 0; d < scanner.size(); d++) {
    if (t_detector.approx_equal_dihedral(scanner[d], epsilon))
      return d;
  }
  return -1;
}

#endif  // FIND_SYMMETRY_H
