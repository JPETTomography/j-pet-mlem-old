#include "util/test.h"

#include "common/types.h"

#include "2d/strip/gausian_kernel.h"
#include "2d/strip/kernel.h"
#include "common/types.h"
#include "2d/geometry/vector.h"

using Kernel = PET2D::Strip::GaussianKernel<F>;
using Vector = PET2D::Vector<F>;

#include "2d/strip/kernel_testing_toolbox.h"

using namespace PET2D::Strip::Testing;

TEST("strip/events") {

  Event<F> evt(0.1, -.2, M_PI / 3);

  CHECK(evt.x == Approx(0.1));
  CHECK(evt.y == Approx(-0.2));
  CHECK(evt.theta == Approx(M_PI / 3));
  CHECK(evt.tan == Approx(tan(M_PI / 3)));

  FrameEvent<F> fevnt(-.2, .1, 0.3);

  CHECK(fevnt.zup == Approx(-0.2));
  CHECK(fevnt.zdn == Approx(0.1));
  CHECK(fevnt.dl == Approx(0.3));
}

TEST("strip/event/conversion/Event-FrameEvent-Event") {
  F R = 0.4;
  Event<F> evt(0.1, -0.2, M_PI / 6);

  Event<F> conv(FrameEvent<F>(evt, R), R);

  CHECK(conv.x == Approx(evt.x));
  CHECK(conv.y == Approx(evt.y));
  CHECK(conv.theta == Approx(evt.theta));
  CHECK(conv.tan == Approx(evt.tan));
}

TEST("strip/event/diff") {

  FrameEvent<F> lhs(.1, .2, 0.3), rhs(-.1, .3, 0.2);

  auto diff = lhs - rhs;

  CHECK(diff.x == Approx(lhs.zup - rhs.zup));
  CHECK(diff.y == Approx(lhs.zdn - rhs.zdn));
  CHECK(diff.z == Approx(lhs.dl - rhs.dl));
}

TEST("strip/vector/diagonalmul") {

  Vector3D<F> diag(0.3, 0.1, 0.2);
  Vector3D<F> vec(2., 2.5, .3);

  auto res = diagonal_product(diag, vec);

  CHECK(res == Approx(2 * 0.3 * 2 + 2.5 * 0.1 * 2.5 + 0.3 * .2 * 0.3));
}

TEST("strip/gauss") {
  Vector3D<F> diag(0.3, 0.1, 0.2);
  Vector3D<F> vec(0.1, 0.2, .3);

  auto res = gauss(diag, vec);
  CHECK(res == Approx(std::exp(-0.5 * (0.3 * 0.1 * 0.1 + 0.1 * 0.2 * 0.2 +
                                       0.2 * 0.3 * 0.3))));
}

TEST("strip/sensitivity") {
  F L = 0.5;
  CHECK(sensitivity(FrameEvent<F>(0, 0, 0), L) == Approx(1.0));
  CHECK(sensitivity(FrameEvent<F>(0, L / 2, 0), L) == Approx(1.0));
  CHECK(sensitivity(FrameEvent<F>(-L / 2, L / 2, 0), L) == Approx(1.0));

  F epsilon = 1e-5;
  CHECK(sensitivity(FrameEvent<F>(0, L / 2 + epsilon, 0), L) == Approx(0.0));
  CHECK(sensitivity(FrameEvent<F>(-L / 2 - epsilon, L / 2, 0), L) ==
        Approx(0.0));
}

TEST("strip/integral") {
  F sz = 0.01;
  F sdl = 0.04;
  F R = 0.4;
  F L = 0.5;

  Vector3D<F> diag(1 / (sz * sz), 1 / (sz * sz), 1 / (sdl * sdl));
  FrameEvent<F> exact(Event<F>(0, 0, 0.2), R);
  F d = 0.01;
  F z_lim = 0.25;
  F dl_lim = 0.8;
  double sum = 0.0;
  for (F zup = -z_lim; zup <= z_lim; zup += d) {
    for (F zdn = -z_lim; zdn <= z_lim; zdn += d) {
      for (F dl = -dl_lim; dl <= dl_lim; dl += d) {
        sum += weight(diag, FrameEvent<F>(zup, zdn, dl), exact, L);
      }
    }
  }
  std::cout << "w integrated : " << sum * d * d * d << "\n";
}

TEST("strip/integral/event") {
  F sz = 0.01;
  F sdl = 0.04;
  F R = 0.43;
  F L = 0.5;

  Vector3D<F> diag(1 / (sz * sz), 1 / (sz * sz), 1 / (sdl * sdl));
  FrameEvent<F> exact(Event<F>(0, 0, 0.2), R);
  F d = 0.01;
  F xy_lim = 0.25;

  double sum = 0.0;
  for (F x = -xy_lim; x < xy_lim; x += d) {
    for (F y = -xy_lim; y < xy_lim; y += d) {
      for (F angle = -M_PI / 4; angle < M_PI / 4; angle += d) {
        sum += weight(diag, FrameEvent<F>(Event<F>(x, y, angle), R), exact, L) *
               4 * R / pow(std::cos(angle), 1.5);
      }
    }
  }
  std::cout << "w integrated (event): " << sum * d * d * d << "\n";
}

TEST("strip/integral/theta") {
  F sz = 0.01;
  F sdl = 0.04;
  F R = 0.43;
  F L = 0.5;

  Vector3D<F> diag(1 / (sz * sz), 1 / (sz * sz), 1 / (sdl * sdl));

  F d = 0.01;
  F z_lim = 0.25;
  F dl_lim = 0.8;
  double sum = 0.0;
  for (F zup = -z_lim; zup <= z_lim; zup += d) {
    for (F zdn = -z_lim; zdn <= z_lim; zdn += d) {
      for (F dl = -dl_lim; dl <= dl_lim; dl += d) {
        sum += theta_integral(diag, FrameEvent<F>(zup, zdn, dl), 0, 0, R, L);
      }
    }
  }
  auto integral = sum * d * d * d;
  auto sens = sensitivity(F(0.0), F(0.0), R, L);
  std::cout << "w integrated theta : " << integral << " /  " << sens << " = "
            << integral / sens << "\n";
}


TEST("strip/gauss_kernel") { Kernel kernel(0.01, 0.04); }

TEST("strip/gauss_kernel/integral") {

  Kernel kernel(0.01, 0.04);

  F y = 0.0;
  F R = 0.43;

  F d = 0.01;
  F dfi = 0.01;

  double sum = 0.0;
  for (F angle = -M_PI / 4; angle < M_PI / 4; angle += dfi) {
    F tan = std::tan(angle);
    F sec = 1 / cos(angle);
    for (F dx = -0.2; dx <= 0.2; dx += d) {
      for (F dy = -0.2; dy <= 0.2; dy += d) {
        sum += kernel(y + dy, tan, sec, R, Vector(dx, dy));
      }
    }
  }
  std::cout << "gauss_kernel integrated : " << sum * d * d * dfi << "\n";
}
