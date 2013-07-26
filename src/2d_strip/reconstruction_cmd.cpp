#include <iostream>
#include <vector>
#include <ctime>
#include <xmmintrin.h>

#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"
#include "util/cmdline_types.h"

#include "phantom.h"
#include "reconstruction.h"
#include "event.h"

using namespace std;

int main(int argc, char* argv[]) {

  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

  try {
    cmdline::parser cl;

#if _OPENMP
    cl.add<int>("n-threads", 't', "number of OpenMP threads", false, 4);
#endif
    cl.add<float>("r-distance", 'r', "R distance between scientilators", false,
                  200.0f);
    cl.add<float>("s-length", 'l', "Scentilator_length", false, 400.0f);
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

    double R_distance = cl.get<float>("r-distance");
    double Scentilator_length = cl.get<float>("s-length");
    int iteration = cl.get<int>("iter");
    double pixel_size = cl.get<float>("p-size");
    int n_pixels = Scentilator_length / pixel_size;
    double sigma = cl.get<float>("s-z");
    double dl = cl.get<float>("s-dl");
    double x = 0.0f;
    double y = 0.0f;
    double a = 50.0f;
    double b = 50.0f;
    double phi = 0.0f;

    phantom<double> test(iteration, n_pixels, pixel_size, R_distance,
                         Scentilator_length, x, y, a, b, phi);

    int n_threads = 4;

    std::clock_t t0, t1;

    t0 = std::clock();

    test.emit_event(n_threads);

    t1 = std::clock();

    std::cout << "Event time: " << (t1 - t0) / 1000 << std::endl;

    std::string fn("test.bin");
    test.save_output(fn);

    std::vector<event<float>> list;

    std::cout << "    REC!!!!    " << std::endl;

    spet_reconstruction<double> reconstruction(R_distance, Scentilator_length,
                                               n_pixels, pixel_size, sigma, dl);
    reconstruction.load_input(fn);

    double x_1 = 0.0;
    double y_1 = 0.0;
    double tan = 1.0;
    double z_u = -200.0f;
    double z_d = 200.0f;
    double t = reconstruction.get_event_tan(z_u, z_d);

    int it = 1;
    reconstruction.reconstruction(it);
    std::cout << "TEST KERNELA DLA y = 0,z = 0, tan = 1,z_u = 200,z_d = -200 " <<std::endl;
    std::cout << "TAN: " << t << std::endl;
    std::pair<int, int> p = reconstruction.in_pixel(x_1, y_1);

    // std::cout << p.first << " " << p.second << std::endl;
    std::cout << "KERNEL: " << reconstruction.kernel(y, tan, p) << std::endl;
  }
  catch (std::string & ex) {
    std::cerr << "error: " << ex << std::endl;
  }
  catch (const char * ex) {
    std::cerr << "error: " << ex << std::endl;
  }

  return 0;
}
