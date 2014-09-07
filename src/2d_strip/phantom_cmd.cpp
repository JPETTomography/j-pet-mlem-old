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
#include "strip_detector.h"

#if _OPENMP
#include <omp.h>
#endif

const double RADIAN = M_PI / 180;

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
    cl.add<int>("n-pixels", 'n', "number of pixels", false);
    cl.add<int>("n-z-pixels", '\0', "number of z pixels", false);
    cl.add<int>("n-y-pixels", '\0', "number of y pixels", false);
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

    double sigma_z = cl.get<double>("s-z");
    double sigma_dl = cl.get<double>("s-dl");
    double emissions = cl.get<double>("emissions");

    // typedef EllipseParameters<double> Ellipse;
    std::vector<PhantomRegion<double>> ellipse_list;

    int n_z_pixels;
    int n_y_pixels;
    if (cl.exist("n-pixels")) {
      n_z_pixels = cl.get<int>("n-pixels");
      n_y_pixels = cl.get<int>("n-pixels");
    } else {
      if (cl.exist("n-z-pixels"))
        n_z_pixels = cl.get<int>("n-z-pixels");
      else
        n_z_pixels = (int)std::ceil(scintillator_length / pixel_size);

      if (cl.exist("n-y-pixels"))
        n_y_pixels = cl.get<int>("nyz-pixels");
      else
        n_y_pixels = (int)std::ceil(2 * R_distance / pixel_size);
    }
    for (auto& fn : cl.rest()) {
      std::ifstream infile(fn);
      std::string line;
      while (std::getline(infile, line)) {
        std::istringstream iss(line);
        double x, y, a, b, angle, acceptance;

        // on error
        if (!(iss >> x >> y >> a >> b >> angle >> acceptance))
          break;

        Ellipse<double> el(x, y, a, b, angle * RADIAN);

        std::cout << el.x() << " " << el.y() << " " << el.a() << " " << el.b()
                  << " " << el.angle() << " " << el.A() << " " << el.B() << " "
                  << el.C() << std::endl;
        PhantomRegion<double> region(el, acceptance);

        ellipse_list.push_back(region);
      }
    }

    Phantom<StripDetector<double>, double> phantom(
        StripDetector<double>(R_distance,
                              scintillator_length,
                              n_y_pixels,
                              n_z_pixels,
                              pixel_size,
                              pixel_size,
                              sigma_z,
                              sigma_dl),
        ellipse_list);

    phantom(emissions);

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
