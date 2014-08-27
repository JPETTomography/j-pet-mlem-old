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
    cl.footer("phantom_description");

    cl.add<cmdline::path>(
        "output", 'o', "output events file", false, "phantom.bin");
    cl.add<double>(
        "r-distance", 'r', "R distance between scientilators", false, 500);
    cl.add<double>("s-length", 'l', "scintillator length", false, 1000);
    cl.add<double>("p-size", 'p', "pixel size", false, 5);
    cl.add<int>("n-pixels", 'n', "number of pixels", false, 200);
    cl.add<double>("s-z", 's', "Sigma z error", false, 10);
    cl.add<double>("s-dl", 'd', "Sigma dl error", false, 63);
    cl.add<double>("emissions", 'e', "number of emissions", false, 500000);
#if _OPENMP
    cl.add<int>("n-threads", 'T', "number of OpenMP threads", false, 4);
#endif

    cl.parse_check(argc, argv);

    if (!cl.rest().size()) {
      throw(
          "at least one input phantom description file expected, "
          "consult --help");
    }

#if _OPENMP
    if (cl.exist("n-threads")) {
      omp_set_num_threads(cl.get<int>("n-threads"));
    }
#endif

    double R_distance = cl.get<double>("r-distance");
    double scintillator_length = cl.get<double>("s-length");
    double pixel_size = cl.get<double>("p-size");
    double n_pixels = scintillator_length / pixel_size;
    double sigma_z = cl.get<double>("s-z");
    double sigma_dl = cl.get<double>("s-dl");
    double emissions = cl.get<double>("emissions");

    typedef EllipseParameters<double> Ellipse;
    std::vector<Ellipse> ellipse_list;

    double normalized_acceptance = 0;

    for (auto& fn : cl.rest()) {
      std::ifstream infile(fn);
      std::string line;
      while (std::getline(infile, line)) {
        std::istringstream iss(line);
        double x, y, a, b, angle, acceptance;

        // on error
        if (!(iss >> x >> y >> a >> b >> angle >> acceptance))
          break;

        Ellipse el(x, y, a, b, angle, acceptance);

        normalized_acceptance += acceptance;

        std::cout << el.x << " " << el.y << " " << el.a << " " << el.b << " "
                  << el.angle << " " << el.n_emissions << std::endl;

        ellipse_list.push_back(el);
      }
    }

    for (auto& el : ellipse_list) {
      el.n_emissions *= emissions / normalized_acceptance;
      std::cout << el.n_emissions << std::endl;
    }

    Phantom<Ellipse::F> phantom(ellipse_list,
                                n_pixels,
                                pixel_size,
                                R_distance,
                                scintillator_length,
                                sigma_z,
                                sigma_dl);

    phantom();

    auto output = cl.get<cmdline::path>("output");
    obstream out(output);
    phantom >> out;

    png_writer png(output.wo_ext() + ".png");
    phantom.output_bitmap(png);

    png_writer png_true(output.wo_ext() + "_true.png");
    phantom.output_bitmap(png_true, true);
  } catch (std::string& ex) {
    std::cerr << "error: " << ex << std::endl;
  } catch (const char* ex) {
    std::cerr << "error: " << ex << std::endl;
  }

  return 0;
}
