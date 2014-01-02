// PET Reconstruction
// Authors:
//   Piotr Bialas    <piotr.bialas@uj.edu.pl>
//   Jakub Kowal     <jakub.kowal@uj.edu.pl>
//   Adam Strzelecki <adam.strzlecki@uj.edu.pl>
//
// Based on:
//   "Implementing and Accelerating the EM Algorithm for Positron Emission
// Tomography"
//   by Linda Kaufman

#ifdef __SSE3__
#include <xmmintrin.h>
#include <pmmintrin.h>
#endif

#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"
#include "util/cmdline_types.h"
#include "reconstruction.h"
#include "util/png_writer.h"

#if _OPENMP
#include <omp.h>
#endif

template <typename Iterator>
void output_vector(std::ostream& out,
                   Iterator start,
                   Iterator stop,
                   int line_length) {
  int index = 0;
  for (; start != stop; ++start, ++index) {
    out << *start << " ";
    if (index % line_length == 0)
      out << "\n";
  }
}

int main(int argc, char* argv[]) {

#ifdef __SSE3__
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
#endif

  try {
    cmdline::parser cl;

#if _OPENMP
    cl.add<int>("n-threads", 't', "number of OpenMP threads", false);
#endif
    cl.add<cmdline::string>("system", 's', "system matrix file", true);
    cl.add<cmdline::string>("mean", 'm', "mean file", true);
    cl.add<int>("iterations", 'n', "number of iterations", false, 0);
    cl.add<int>("i-blocks", 'i', "number of iteration blocks", false, 1);
    cl.add<cmdline::string>("output", 'o', "output reconstruction", false);
    cl.add<double>("threshold", 0, "discretisation treshold", false, 0.0);

    // additional options
    cl.add("verbose", 'v', "prints the iterations information on std::out");

    cl.parse_check(argc, argv);

#if _OPENMP
    if (cl.exist("n-threads")) {
      omp_set_num_threads(cl.get<int>("n-threads"));
    }
#endif

    auto verbose = cl.exist("verbose");

    ibstream in_matrix(cl.get<cmdline::string>("system"));
    if (!in_matrix.is_open())
      throw("cannot open input file: " + cl.get<cmdline::string>("system"));
    Reconstruction<>::Matrix matrix(in_matrix);
    matrix = matrix.to_full();

    std::ifstream in_means(cl.get<cmdline::string>("mean"));
    if (!in_means.is_open())
      throw("cannot open input file: " + cl.get<cmdline::string>("mean"));

    int n_i_blocks = cl.get<int>("i-blocks");
    Reconstruction<> reconstruction(cl.get<int>("iterations"),
                                    matrix,
                                    in_means,
                                    cl.get<double>("threshold"));

    if (verbose) {
      std::cerr << "reconstruction:" << std::endl;
#if _OPENMP
      std::cerr << " threads       = " << omp_get_max_threads() << std::endl;
#endif
      std::cerr << " pixels in row = " << matrix.n_pixels_in_row() << std::endl;
      std::cerr << " TOF positions = " << matrix.n_tof_positions() << std::endl;
      std::cerr << " emissions     = " << matrix.n_emissions() << std::endl;
      std::cerr << " detectors     = " << matrix.n_detectors() << std::endl;
    }

    auto n_pixels_in_row = reconstruction.n_pixels_in_row();
    auto total_n_pixels = n_pixels_in_row * n_pixels_in_row;
    Reconstruction<>::Output rho(total_n_pixels, 0.0);
    Reconstruction<>::Output rho_detected(total_n_pixels, 0.0);

    // no output, just make reconstruction in place and exit!
    if (!cl.exist("output")) {
      for (int i = 0; i < n_i_blocks; ++i) {
        reconstruction.emt(cl.get<int>("iterations"));
        rho = reconstruction.rho();
      }
      return 0;
    }

    std::string fn;
    std::string fn_ext;
    std::string fn_wo_ext;

    fn = cl.get<cmdline::string>("output");
    auto it_fn_sep = fn.find_last_of("\\/");
    auto it_fn_ext = fn.find_last_of(".");
    fn_wo_ext = fn.substr(
        0,
        it_fn_ext != std::string::npos &&
                (it_fn_sep == std::string::npos || it_fn_sep < it_fn_ext)
            ? it_fn_ext
            : std::string::npos);
    fn_ext = fn.substr(it_fn_ext != std::string::npos ? it_fn_ext : fn.size(),
                       fn.size());

    std::ofstream out;
    std::ofstream out_detected;

    out.open(fn);
    out_detected.open(fn_wo_ext + "_detected" + fn_ext);

    for (int i = 0; i < n_i_blocks; ++i) {
      reconstruction.emt(cl.get<int>("iterations"));
      rho = reconstruction.rho();
      rho_detected = reconstruction.rho_detected();
      output_vector(out, rho.begin(), rho.end(), n_pixels_in_row);
      output_vector(out_detected,
                    rho_detected.begin(),
                    rho_detected.end(),
                    n_pixels_in_row);
    }

    out.close();
    out_detected.close();

    // output reconstruction PNG
    png_writer png(fn_wo_ext + ".png");
    png.write_header<>(n_pixels_in_row, n_pixels_in_row);

    double output_max = 0.0;
    for (auto it = rho.begin(); it != rho.end(); ++it) {
      output_max = std::max(output_max, *it);
    }

    auto output_gain =
        static_cast<double>(std::numeric_limits<uint8_t>::max()) / output_max;

    for (int y = n_pixels_in_row - 1; y >= 0; --y) {
      uint8_t row[n_pixels_in_row];
      for (auto x = 0; x < n_pixels_in_row; ++x) {
        row[x] = std::numeric_limits<uint8_t>::max() -
                 output_gain * rho[y * n_pixels_in_row + x];
      }
      png.write_row(row);
    }

    return 0;
  }
  catch (std::string& ex) {
    std::cerr << "error: " << ex << std::endl;
  }
  catch (const char* ex) {
    std::cerr << "error: " << ex << std::endl;
  }
  return 1;
}
