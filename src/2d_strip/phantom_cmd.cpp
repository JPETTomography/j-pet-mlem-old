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
#include "options.h"

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

    add_phantom_options(cl);

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

    double emissions = cl.get<double>("emissions");

    std::vector<PhantomRegion<double>> ellipse_list;

    StripDetector<double> detector = strip_detector_from_options<double>(cl);

    std::cerr << detector.n_z_pixels << "x" << detector.n_y_pixels << std::endl;
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

        std::cout << el.x << " " << el.y << " " << el.a << " " << el.b << " "
                  << el.angle << " " << el.A << " " << el.B << " " << el.C
                  << std::endl;
        PhantomRegion<double> region(el, acceptance);

        ellipse_list.push_back(region);
      }
    }

    Phantom<StripDetector<double>, double> phantom(detector, ellipse_list);
    std::cerr << "detector " << detector.size_y << " " << detector.tl_y_half_h
              << std::endl;

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
