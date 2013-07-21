#include <iostream>
#include <vector>
#include <ctime>
#include <xmmintrin.h>

#include "phantom.h"
#include "reconstruction.h"
#include "event.h"


using namespace std;

int main() {

  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

  float R_distance = 400.0f;
  float Scentilator_length = 400.0f;
  int iteration = 100000;
  float pixel_size = 5.0f;
  int n_pixels = Scentilator_length / pixel_size;
  float x = 0.0f;
  float y = 0.0f;
  float a = 10.0f;
  float b = 63.0f;
  float phi = 0.0f;

  phantom<> test(iteration,
                 n_pixels,
                 pixel_size,
                 R_distance,
                 Scentilator_length,
                 x,
                 y,
                 a,
                 b,
                 phi);

  int n_threads = 4;

  std::clock_t t0, t1;

  t0 = std::clock();

  test.emit_event(n_threads);

  t1 = std::clock();

  std::cout << "Event time: " << (t1 - t0) / 1000 << std::endl;

  std::string fn("test.bin");
  test.save_output(fn);
  test.load_input(fn);

  float sigma = 10.0f;
  float dl = 63.0f;
  float gamma = 0.0f;
  std::vector<event<float>> list;

  std::cout << "    REC!!!!    " << std::endl;

  spet_reconstruction<> reconstruction(
      R_distance, Scentilator_length, n_pixels, pixel_size, sigma, dl, gamma);
  // reconstruction.load_input(fn);

  float x_1 = 0.0;
  float y_1 = 0.0;
  float tan = 1.0;
  std::pair<int, int> p = reconstruction.in_pixel(x_1, y_1);

  std::cout << p.first << " " << p.second << std::endl;
  std::cout << "KERNEL: " << reconstruction.kernel(y, tan, p) << std::endl;

  return 0;
}
