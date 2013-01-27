// PET Reconstruction
// Authors:
//   Piotr Bialas    <piotr.bialas@uj.edu.pl>
//   Jakub Kowal     <jakub.kowal@uj.edu.pl>
//   Adam Strzelecki <adam.strzlecki@uj.edu.pl>
//
// Based on:
//   "Implementing and Accelerating the EM Algorithm for Positron Emission Tomography"
//   by Linda Kaufman

#include <xmmintrin.h>
#include <pmmintrin.h>

#include <cmdline.h>
#include "cmdline_types.h"
#include "bstream.h"
#include "svg_ostream.h"
#include "cmdline_types.h"
#include "reconstruction.h"
#include "png_writer.h"

#if _OPENMP
#include <omp.h>
#endif

template <typename Iterator>
void write_text_from_vector(
    std::ostream& out, Iterator start, Iterator stop, int line_length) {
  int index = 0;
  for (; start != stop; ++start, ++index) {
    out << *start << " ";
    if (index % line_length == 0)
      out << "\n";
  }
}

int main(int argc, char* argv[]) {

  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);

  try {
    cmdline::parser cl;

#if _OPENMP
    cl.add<size_t>("n-threads", 't', "number of OpenMP threads", false);
#endif
    cl.add<cmdline::string>("system", 's', "system matrix file", true);
    cl.add<cmdline::string>("mean", 'm', "mean file", true);
    cl.add<size_t>("iterations", 'n', "number of iterations", false, 0);
    cl.add<int>("i-blocks", 'i', "number of iteration blocks", false, 1);
    cl.add<cmdline::string>("output", 'o', "output reconstruction", false);
    cl.add<double>("threshold", '\0', "discretisation treshold", false, 0.0);

    cl.parse_check(argc, argv);
    std::string fn;
    std::string fn_ext;
    std::string fn_wo_ext;
    bool do_output = cl.exist("output");

    if (do_output) {

      fn = cl.get<cmdline::string>("output");
      auto it_fn_sep = fn.find_last_of("\\/");
      auto it_fn_ext = fn.find_last_of(".");
      fn_wo_ext =
          fn.substr(0, it_fn_ext != std::string::npos &&
                    (it_fn_sep == std::string::npos || it_fn_sep < it_fn_ext) ?
                        it_fn_ext : std::string::npos);
      fn_ext = fn.substr(it_fn_ext != std::string::npos ? it_fn_ext : fn.size(),
                         fn.size());
    }

#if _OPENMP
    if (cl.exist("n-threads")) {
      omp_set_num_threads(cl.get<size_t>("n-threads"));
    }
#endif

    int n_i_blocks = cl.get<int>("i-blocks");
    reconstruction<> reconstructor(cl.get<size_t>("iterations"),
                                   cl.get<cmdline::string>("system"),
                                   cl.get<cmdline::string>("mean"),
                                   cl.get<double>("threshold"));

    auto n_pixels = reconstructor.get_n_pixels();

    auto total_n_pixels = n_pixels * n_pixels;
    reconstruction<>::output_type rho(total_n_pixels, 0.0);
    reconstruction<>::output_type rho_detected(total_n_pixels, 0.0);

    std::ofstream out;
    std::ofstream out_detected;
    if (do_output) {
      out.open(fn);
      out_detected.open(fn_wo_ext + "_detected" + fn_ext);
    }
    for (int i = 0; i < n_i_blocks; ++i) {
      reconstructor.emt(cl.get<size_t>("iterations"));
      rho = reconstructor.rho();
      rho_detected = reconstructor.rho_detected();
      if (do_output) {
        write_text_from_vector(out, rho.begin(), rho.end(), n_pixels);
        write_text_from_vector(
            out_detected, rho_detected.begin(), rho_detected.end(), n_pixels);
      }
    }
    out.close();
    out_detected.close();

    if (!do_output)
      return 0;

    // output reconstruction PNG
    png_writer png(fn_wo_ext + ".png");
    png.write_header<>(n_pixels, n_pixels);

    double output_max = 0.0;
    for (auto it = rho.begin(); it != rho.end(); ++it) {
      output_max = std::max(output_max, *it);
    }

    auto output_gain =
        static_cast<double>(std::numeric_limits<uint8_t>::max()) / output_max;

    for (int y = n_pixels - 1; y >= 0; --y) {
      uint8_t row[n_pixels];
      for (auto x = 0; x < n_pixels; ++x) {
        row[x] = std::numeric_limits<uint8_t>::max() -
                 output_gain * rho[LOCATION(x, y, n_pixels)];
      }
      png.write_row(row);
    }

    return 0;

  }
  catch (std::string & ex) {
    std::cerr << "error: " << ex << std::endl;
  }
  catch (const char * ex) {
    std::cerr << "error: " << ex << std::endl;
  }

}
