#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"

#include "2d/strip/gaussian_kernel.h"

#include "common/types.h"
#include "2d/geometry/vector.h"

#include "2d/strip/kernel_testing_toolbox.h"

using Kernel = PET2D::Strip::GaussianKernel<F>;
using Vector = PET2D::Vector<F>;
using namespace PET2D::Strip;
using namespace PET2D::Strip::Testing;

void strip_integral() {
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

void strip_integral_theta() {
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

void strip_integral_event() {
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
    for (F y = -.4; y < 0.4; y += d) {
      for (F angle = -M_PI / 4; angle < M_PI / 4; angle += d) {
        sum += weight(diag, FrameEvent<F>(Event<F>(x, y, angle), R), exact, L) *
               4 * R / pow(std::cos(angle), 3);
      }
    }
  }
  std::cout << "w integrated (event): " << sum * d * d * d << "\n";
}

void strip_integral_theta_event() {

  F sz = 0.01;
  F sdl = 0.04;
  F R = 0.43;
  F L = 0.5;

  Vector3D<F> diag(1 / (sz * sz), 1 / (sz * sz), 1 / (sdl * sdl));

  F d = 0.01;
  F xy_lim = 0.25;

  double sum = 0.0;
  for (F x = -xy_lim; x < xy_lim; x += d) {
    for (F y = -0.4; y < 0.4; y += d) {
      for (F angle = -M_PI / 4; angle < M_PI / 4; angle += d) {
        sum += theta_integral(
                   diag, FrameEvent<F>(Event<F>(x, y, angle), R), 0, 0, R, L) *
               4 * R / pow(std::cos(angle), 3);
      }
    }
  }

  auto integral = sum * d * d * d;
  auto sens = sensitivity(F(0.0), F(0.0), R, L);
  std::cout << "w integrated theta (event): " << integral << " /  " << sens
            << " = " << integral / sens << "\n";
}

void strip_gauss_kernel() {

  F R = 0.43;
  F L = 0.5;
  F sz = 0.01;
  F sdl = 0.04;

  Kernel kernel(sz, sdl);

  F tx = 0.01;
  F ty = 0.03;
  F tangle = 0.0;

  F x = 0.005;
  F y = 0.0;

  Vector3D<F> diag(1 / (sz * sz), 1 / (sz * sz), 1 / (sdl * sdl));

  std::cout << " theta = "
            << theta_integral(
                   diag, FrameEvent<F>(Event<F>(tx, ty, tangle), R), x, y, R, L)
            << std::endl;
  std::cout << "kernel = " << kernel(ty,
                                     std::tan(tangle),
                                     1 / std::cos(tangle),
                                     R,
                                     Vector(x - tx, y - ty))
            << std::endl;
}

void strip_gauss_kernel_integral() {
  F L = 0.5;
  Kernel kernel(0.01, 0.04);

  F x = 0.0;
  F y = 0.0;
  F R = 0.43;

  F d = 0.002;
  F dfi = 0.01;

  double sum = 0.0;

  for (F tx = -0.25; tx <= 0.25; tx += d) {
    for (F ty = -0.4; ty <= 0.4; ty += d) {
      auto theta_lim = theta_min_max(tx, ty, R, L);
      for (F tangle = theta_lim.first; tangle < theta_lim.second;
           tangle += dfi) {
        F tan = std::tan(tangle);
        F sec = 1 / cos(tangle);
        sum += kernel(ty, tan, sec, R, Vector(x - tx, y - ty)) * 4 * R /
               pow(std::cos(tangle), 3.0);
      }
    }
  }
  auto sens = sensitivity(F(0.0), F(0.0), R, L);
  auto integral = sum * d * d * dfi;

  std::cout << "gauss_kernel integrated : " << integral << " / " << sens
            << " = " << integral / sens << "\n";
}

int main(int argc, char* argv[]) {

  cmdline::parser cl;
  cl.parse_check(argc, argv);

  strip_integral();
  strip_integral_theta();
  strip_integral_event();
  strip_integral_theta_event();
  strip_gauss_kernel();
  strip_gauss_kernel_integral();
}
