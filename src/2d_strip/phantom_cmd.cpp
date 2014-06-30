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

#include "phantom.h"

#if _OPENMP
#include <omp.h>
#endif

int main(int argc, char* argv[]) {

#if SSE_FLUSH
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif

  try {
    cmdline::parser cl;

#if _OPENMP
    cl.add<int>("n-threads", 't', "number of OpenMP threads", false, 4);
#endif
    cl.add<cmdline::string>("output", 'o', "events file", false, "phantom.bin");
    cl.add<double>(
        "r-distance", 'r', "R distance between scientilators", false, 500);
    cl.add<double>("s-length", 'l', "Scentilator_length", false, 1000);
    cl.add<double>("p-size", 'p', "Pixel size", false, 5);
    cl.add<int>("n-pixels", 'n', "Number of pixels", false, 200);
    cl.add<int>("iter", 'i', "number of iterations", false, 1);
    cl.add<double>("s-z", 's', "Sigma z error", false, 10);
    cl.add<double>("s-dl", 'd', "Sigma dl error", false, 63);
    cl.add<double>("gm", 'g', "Gamma error", false, 0);
    cl.add<double>("emmisions", 'e', "number of emissions", false, 500000);

    cl.parse_check(argc, argv);

#if _OPENMP
    if (cl.exist("n-threads")) {
      omp_set_num_threads(cl.get<int>("n-threads"));
    }
#endif

    double R_distance = cl.get<double>("r-distance");
    double Scentilator_length = cl.get<double>("s-length");
    double pixel_size = cl.get<double>("p-size");
    double n_pixels = Scentilator_length / pixel_size;
    double sigma = cl.get<double>("s-z");
    double dl = cl.get<double>("s-dl");
    double emmisions = cl.get<double>("emmisions");

    typedef EllipseParameters<double> Ellipse;
    std::vector<Ellipse> ellipse_list;

    double normalized_acc = 0;

    for (auto& fn : cl.rest()) {
      std::ifstream infile(fn);
      std::string line;
      while (std::getline(infile, line)) {
        std::istringstream iss(line);
        double x, y, a, b, angle, acc;

        // on error
        if (!(iss >> x >> y >> a >> b >> angle >> acc))
          break;

        Ellipse el(x, y, a, b, angle, acc);

        normalized_acc += acc;

        std::cout << el.x << " " << el.y << " " << el.a << " " << el.b << " "
                  << el.angle << " " << el.iter << std::endl;

        ellipse_list.push_back(el);
      }
    }

    for (auto& e : ellipse_list) {
      e.iter = e.iter * (emmisions / normalized_acc);
      std::cout << e.iter << std::endl;
    }

    Phantom<Ellipse::F> test(ellipse_list,
                             n_pixels,
                             pixel_size,
                             R_distance,
                             Scentilator_length,
                             sigma,
                             dl);

    test.emit_event();

    obstream out(cl.get<cmdline::string>("output"));
    test >> out;
  } catch (std::string& ex) {
    std::cerr << "error: " << ex << std::endl;
  } catch (const char* ex) {
    std::cerr << "error: " << ex << std::endl;
  }

  return 0;
}
