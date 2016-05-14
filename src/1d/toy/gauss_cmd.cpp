#include <random>
#include <vector>
#include <iostream>

#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"

#include "common/types.h"

#include "1d/toy/gauss_kernel.h"

F s(F x) {
  if (x <= -2)
    return 0;
  if (x <= -1)
    return x + 2;
  if (x <= 1)
    return 1;
  if (x <= 2)
    return -x + 2;
  return 0;
}

using Kernel = PET1D::Toy::GaussKernel<F>;

int main(int argc, char* argv[]) {

  cmdline::parser cl;
  cl.add<int>("emissions", 'e', "number of emissions", false, 0);
  cl.add<F>("dx", '\0', "delta x", false, 0.01);
  cl.parse_check(argc, argv);

  size_t n = cl.get<int>("emissions");
  F dx = cl.get<F>("dx");
  int n_x = 4 / dx;

  F sigma = 0.25;
  Kernel kernel(sigma);

  std::mt19937_64 rng;
  std::uniform_real_distribution<F> uni;
  std::uniform_real_distribution<F> unix(-2, 2);
  std::normal_distribution<F> error(0, sigma);

  std::vector<F> xs, ys;
  size_t n_samples = 0;
  for (size_t i = 0; i < n; ++i) {
    F x = unix(rng);
    F r = uni(rng);
    if (r < s(x)) {
      n_samples++;
      xs.push_back(x);
      ys.push_back(x + error(rng));
    }
  }
  std::cout << "# " << n_samples << "\n";

  std::vector<F> prob_ys(n_samples);
#if _OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < (int)n_samples; ++i) {
    F p = 0;
    for (F x = -2; x <= 2.0; x += dx) {
      p += kernel(ys[i], x) * s(x);
    }
    prob_ys[i] = p * dx;
  }

  std::vector<F> xv(n_x);
  if (n > 0) {
#if _OPENMP
#pragma omp parallel for schedule(dynamic, 4)
#endif
    for (int j = 0; j < n_x; ++j) {
      F x = -2 + j * dx;
      F sum = 0.0;
      for (size_t i = 0; i < n_samples; ++i) {
        sum += kernel(ys[i], x) / prob_ys[i];
      }
      xv[j] = sum;
    }
    for (int j = 0; j < n_x; ++j) {
      F x = -2 + j * dx;
      std::cout << x << " " << xv[j] * dx << "\n";
    }
  }
}
