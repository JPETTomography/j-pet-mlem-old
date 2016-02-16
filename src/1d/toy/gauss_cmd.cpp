#include <random>
#include <vector>
#include <iostream>

#include "common/types.h"

#include "1d/toy/gauss_kernel.h"

F s(F x) {
  if (x >= 2)
    return 0;
  if (x >= 0)
    return -0.5 * x + 1;
  if (x >= -2)
    return 0.5 * x + 1;
  return 0;
}

using Kernel = PET1D::Toy::GaussKernel<F>;

int main(int argc, const char* argv[]) {

  F sigma = 0.25;
  Kernel kernel(sigma);

  std::mt19937_64 rng;
  std::uniform_real_distribution<F> uni;
  std::uniform_real_distribution<F> unix(-2, 2);
  std::normal_distribution<F> error(0, sigma);

  size_t n = 16000000;
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

  std::vector<F> prob_ys;

  F dx = 0.002;
  size_t count = 0;
  for (auto y : ys) {
    F p = 0;
    for (F x = -2; x <= 2.0; x += dx) {
      p += kernel(y, x) * s(x);
    }
    prob_ys.push_back(p * dx);
    count++;
    if (count % 100000 == 0)
      std::cerr << count << "\n";
  }

  for (F x = -2.0; x < 2.0; x += dx) {
    F sum = 0.0;
    for (size_t i = 0; i < n_samples; ++i) {
      sum += kernel(ys[i], x) / prob_ys[i];
    }

    std::cout << x << " " << sum * dx << "\n";
  }
}
