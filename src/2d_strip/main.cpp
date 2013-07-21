#include <iostream>
#include <vector>
#include <ctime>

#include "phantom.h"
#include "spet_reconstruction.h"
#include "data_structures.h"
#include "omp.h"

using namespace std;

int main() {
  float R_distance = 10.0f;
  float Scentilator_length = 10.0f;
  int iteration = 1000000;
  float pixel_size = 0.5f;
  int n_pixels = 40;
  float x = 0.0f;
  float y = 0.0f;
  float a = 4.0f;
  float b = 3.0f;
  float phi = 0.0f;

  phantom<> test(iteration, n_pixels, pixel_size, R_distance,
                 Scentilator_length, x, y, a, b, phi);

  int n_threads = 4;

  std::clock_t t0, t1;

  t0 = std::clock();

  test.emit_event(n_threads);

  t1 = std::clock();

  std::cout << "Event time: " << (t1 - t0) / 1000 << std::endl;

  std::string fn("output.bin");
  test.save_output(fn);
  test.load_input(fn);

  float sigma = 1;
  float dl = 1;
  float gamma = 0;
  std::vector<event<float> > list;

  spet_reconstruction<> reconstruction(R_distance, Scentilator_length, n_pixels,
                                       pixel_size, sigma, dl, gamma);
  reconstruction.load_input(fn);

  float x_1 = -2.2;
  float y_1 = -2.2;
  float tan = 0.1;
  std::pair<int, int> p = reconstruction.in_pixel(x_1, y_1);
  reconstruction.kernel(y, tan, p);

  std::cout << p.first << " " << p.second << std::endl;

  return 0;
}
