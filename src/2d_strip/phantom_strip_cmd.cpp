#include <iostream>
#include <vector>
#include <ctime>
#include <fstream>
#include <sstream>
#include <string>

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

#if _OPENMP
#include <omp.h>
#endif

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

    float R_distance = cl.get<float>("r-distance");
    float Scentilator_length = cl.get<float>("s-length");
    float pixel_size = cl.get<float>("p-size");
    float n_pixels = Scentilator_length / pixel_size;
    float sigma = cl.get<float>("s-z");
    float dl = cl.get<float>("s-dl");
    float emmisions = cl.get<float>("emmisions");

    std::ifstream infile("phantom");

    std::vector<ellipse_parameters<float>> ellipse_list;
    ellipse_parameters<float> el;

    std::string line;
    while (std::getline(infile, line)) {
      std::istringstream iss(line);
      float x, y, a, b, angle, acc;
      if (!(iss >> x >> y >> a >> b >> angle >> acc)) {
        break;
      }  // error

      el.x = x;
      el.y = y;
      el.a = a;
      el.b = b;
      el.angle = angle;
      el.iter = int(acc * emmisions);

      ellipse_list.push_back(el);
    }

    Phantom<float> test(ellipse_list,
                        n_pixels,
                        pixel_size,
                        R_distance,
                        Scentilator_length,
                        sigma,
                        dl);

    test.emit_event();

    obstream out("phantom.bin");
    test >> out;
  }
  catch (std::string& ex) {
    std::cerr << "error: " << ex << std::endl;
  }
  catch (const char* ex) {
    std::cerr << "error: " << ex << std::endl;
  }

  return 0;
}