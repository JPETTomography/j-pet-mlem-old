#include <iostream>
#include <vector>
#include <ctime>

#if SSE_FLUSH
#include <xmmintrin.h>
#endif

#include "cmdline.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"
#include "util/cmdline_types.h"
#include "util/png_writer.h"

#include "flags.h"
#include "phantom.h"
#include "reconstruction.h"
#include "event.h"

using namespace std;

int main(int argc, char* argv[]) {

#if SSE_FLUSH
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif

  try {
    cmdline::parser cl;

#if _OPENMP
    cl.add<int>("n-threads", 't', "number of OpenMP threads", false, 4);
#endif
    cl.add<float>(
        "r-distance", 'r', "R distance between scientilators", false, 500.0f);
    cl.add<float>("s-length", 'l', "Scentilator_length", false, 1000.0f);
    cl.add<float>("p-size", 'p', "Pixel size", false, 5.0f);
    cl.add<int>("n-pixels", 'n', "Number of pixels", false, 200);
    cl.add<int>("iter", 'i', "number of iterations", false, 1);
    cl.add<float>("s-z", 's', "Sigma z error", false, 10.0f);
    cl.add<float>("s-dl", 'd', "Sigma dl error", false, 63.0f);
    cl.add<float>("gm", 'g', "Gamma error", false, 0.f);
    cl.add<float>("emmisions", 'e', "number of emissions", false, 500000);

    cl.parse_check(argc, argv);

#if _OPENMP
    if (cl.exist("n-threads")) {
      omp_set_num_threads(cl.get<int>("n-threads"));
    }
#endif

    double R_distance = cl.get<float>("r-distance");
    double Scentilator_length = cl.get<float>("s-length");
    double pixel_size = cl.get<float>("p-size");
    double n_pixels = Scentilator_length / pixel_size;
    double sigma = cl.get<float>("s-z");
    double dl = cl.get<float>("s-dl");
    double emmisions = cl.get<float>("emmisions");

    std::vector<ellipse_parameters<double>> ellipse_list;
    ellipse_parameters<double> el;

    el.x = 0.0f;
    el.y = 0.0f;
    el.a = 120.0f;
    el.b = 240.0f;
    el.angle = 0.0;
    el.iter = emmisions;

    ellipse_list.push_back(el);

    el.x = 50.0f;
    el.y = -50.0f;
    el.a = 30.0f;
    el.b = 50.0f;
    el.angle = 340.0;
    el.iter = 200000;

    ellipse_list.push_back(el);
/*
    el.x = -30.0f;
    el.y = -30.0f;
    el.a = 30.0f;
    el.b = 50.0f;
    el.angle = 45.0;
    el.iter = emmisions;

    ellipse_list.push_back(el);

    el.x = 0.0f;
    el.y = 40.0f;
    el.a = 30.0f;
    el.b = 50.0f;
    el.angle = 0.f;
    el.iter = 400000;

    ellipse_list.push_back(el);
*/
    Phantom<double> test(ellipse_list,
                         n_pixels,
                         pixel_size,
                         R_distance,
                         Scentilator_length,
                         sigma,
                         dl);
    int n_threads = cl.get<int>("n-threads");

    test.emit_event(n_threads);

    obstream out("test.bin");
    test >> out;

    Reconstruction<double> reconstruction(cl.get<int>("iter"),
                                          R_distance,
                                          Scentilator_length,
                                          n_pixels,
                                          pixel_size,
                                          sigma,
                                          dl);
    ibstream in("test.bin");

    in >> reconstruction;

    clock_t begin = clock();

    reconstruction();

    clock_t end = clock();

    std::cout << "Time:" << double(end - begin) / CLOCKS_PER_SEC / 4
              << std::endl;
  }
  catch (std::string& ex) {
    std::cerr << "error: " << ex << std::endl;
  }
  catch (const char* ex) {
    std::cerr << "error: " << ex << std::endl;
  }

  return 0;
}
