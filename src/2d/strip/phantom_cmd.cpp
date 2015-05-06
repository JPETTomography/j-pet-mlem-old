/// \page cmd_2d_strip_phantom 2d_strip_phantom
/// \brief 2D Strip PET phantom tool
///
/// Simulates scanner response for given virtual phantom and produces mean file
/// for \ref cmd_2d_strip_reconstruction.
///
/// Example phantom descriptions
/// ----------------------------
/// - Shepp like phantom
///
///   \verbinclude phantom/s_shepp
///
/// - Small Shepp like phantom
///
///   \verbinclude phantom/s_shepp_small
///
/// Authors
/// -------
/// - Adam Strzelecki <adam.strzelecki@uj.edu.pl>
/// - Jakub Kowal     <jakub.kowal@uj.edu.pl>
///
/// Usage
/// -----
/// \verbinclude src/2d/strip/phantom_cmd.txt
///
/// \sa \ref cmd_2d_strip_reconstruction

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
#include "scanner.h"

#if _OPENMP
#include <omp.h>
#endif

using namespace PET2D;
using namespace PET2D::Strip;

const double RADIAN = M_PI / 180;

int main(int argc, char* argv[]) {

#if SSE_FLUSH
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif

  try {
    cmdline::parser cl;
    add_phantom_options(cl);
    cl.parse_check(argc, argv);
    calculate_scanner_options(cl);

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

    auto emissions = cl.get<int>("emissions");
    auto verbose = cl.exist("verbose");

    std::vector<PhantomRegion<double>> ellipse_list;

    Scanner<double, short> scanner(PET2D_STRIP_SCANNER_CL(cl));

    if (verbose) {
      std::cerr << "size: " << scanner.n_z_pixels << "x" << scanner.n_y_pixels
                << std::endl;
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

        if (verbose) {
          std::cout << "ellipse: " << el.center.x << " " << el.center.y << " "
                    << el.a << " " << el.b << " " << el.angle << " " << el.A
                    << " " << el.B << " " << el.C << std::endl;
        }

        PhantomRegion<double> region(el, acceptance);
        ellipse_list.push_back(region);
      }
    }

    Phantom<Scanner<double, short>> phantom(scanner, ellipse_list);

    if (verbose) {
      std::cerr << "scanner: " << scanner.size_y << " " << scanner.tl_y_half_h
                << std::endl;
    }

    phantom(emissions);

    if (verbose) {
      std::cerr << "detected: " << phantom.n_events() << " events" << std::endl;
    }

    auto output = cl.get<cmdline::path>("output");
    auto output_base_name = output.wo_ext();
    util::obstream out(output);
    phantom >> out;

    std::ofstream cfg(output_base_name + ".cfg");
    cfg << cl;

    util::png_writer png(output_base_name + ".png");
    phantom.output_bitmap(png);
    util::obstream bin(output_base_name + "_detected.bin");
    phantom.output_binary(bin, false);

    util::png_writer png_true(output_base_name + "_true.png");
    phantom.output_bitmap(png_true, true);
    util::obstream bin_true(output_base_name + "_true.bin");
    phantom.output_binary(bin_true, true);

  } catch (std::string& ex) {
    std::cerr << "error: " << ex << std::endl;
  } catch (const char* ex) {
    std::cerr << "error: " << ex << std::endl;
  }

  return 0;
}
