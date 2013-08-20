#include <iostream>
#include <vector>
#include <ctime>
#include <xmmintrin.h>

#include "cmdline.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"
#include "util/cmdline_types.h"
#include "util/png_writer.h"

#include "phantom.h"
#include "reconstruction.h"
#include "event.h"

using namespace std;

int main(int argc, char *argv[]) {

  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

  try {
    cmdline::parser cl;

#if _OPENMP
    cl.add<int>("n-threads", 't', "number of OpenMP threads", false, 4);
#endif
    cl.add<float>("r-distance", 'r', "R distance between scientilators", false,
                  500.0f);
    cl.add<float>("s-length", 'l', "Scentilator_length", false, 1000.0f);
    cl.add<float>("p-size", 'p', "Pixel size", false, 5.0f);
    cl.add<int>("n-pixels", 'n', "Number of pixels", false, 80);
    cl.add<int>("iter", 'i', "Reconstruction iterations", false, 20);
    cl.add<float>("s-z", 's', "Sigma z error", false, 10.0f);
    cl.add<float>("s-dl", 'd', "Sigma dl error", false, 63.0f);
    cl.add<float>("gm", 'g', "Gamma error", false, 0.f);

    cl.parse_check(argc, argv);

#if _OPENMP
    if (cl.exist("n-threads")) {
      omp_set_num_threads(cl.get<int>("n-threads"));
    }
#endif

    float R_distance = cl.get<float>("r-distance");
    float Scentilator_length = cl.get<float>("s-length");
    int iteration = 500000;
    float pixel_size = cl.get<float>("p-size");
    int n_pixels = Scentilator_length / pixel_size;
    float sigma = cl.get<float>("s-z");
    float dl = cl.get<float>("s-dl");
    float x = 0.0f;
    float y = 0.0f;
    float a = 30.0f;
    float b = 100.0f;
    float phi = 30.0f;

    Phantom<float> test(iteration, n_pixels, pixel_size, R_distance,
                         Scentilator_length, x, y, a, b, phi);

    int n_threads = 4;

    std::clock_t t0, t1;

    t0 = std::clock();

    test.emit_event(n_threads);

    t1 = std::clock();

    std::cout << "Event time: " << (t1 - t0) / 1000 << std::endl;

    obstream out("test.bin");
    test >> out;

    Reconstruction<float> reconstruction(1, R_distance, Scentilator_length,
                                          n_pixels, pixel_size, sigma, dl);
    ibstream in("test.bin");
    in >> reconstruction;
    //reconstruction.test();
    reconstruction();
  }
  catch (std::string & ex) {
    std::cerr << "error: " << ex << std::endl;
  }
  catch (const char * ex) {
    std::cerr << "error: " << ex << std::endl;
  }

  return 0;
}
