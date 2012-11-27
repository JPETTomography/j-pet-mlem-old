// PET Reconstruction
// Authors:
//   Piotr Bialas    <piotr.bialas@uj.edu.pl>
//   Jakub Kowal     <jakub.kowal@uj.edu.pl>
//   Adam Strzelecki <adam.strzlecki@uj.edu.pl>
//
// Based on:
//   "Implementing and Accelerating the EM Algorithm for Positron Emission Tomography"
//   by Linda Kaufman

#include <cmdline.h>
#include "cmdline_types.h"
#include "bstream.h"
#include "svg_ostream.h"
#include "cmdline_types.h"
#include "pet_rec.h"
#include "png_writer.h"

#if _OPENMP
#include <omp.h>
#endif

int main(int argc, char *argv[]) {

try {
  cmdline::parser cl;

#if _OPENMP
  cl.add<size_t>         ("n-threads", 't', "number of OpenMP threads", false);
#endif
  cl.add<cmdline::string>("system",    's', "system matrix file",       true);
  cl.add<cmdline::string>("mean",      'm', "mean file",                true);
  cl.add<size_t>         ("iterations",'n', "number of iterations",     false, 0);
  cl.add<cmdline::string>("output",    'o', "output reconstruction",    false);

  cl.parse_check(argc, argv);

#if _OPENMP
  if (cl.exist("n-threads")) {
    omp_set_num_threads(cl.get<size_t>("n-threads"));
  }
#endif

  reconstruction<> r(
    cl.get<size_t>("iterations"),
    cl.get<cmdline::string>("system"),
    cl.get<cmdline::string>("mean"));

  auto n_pixels = r.get_n_pixels();
  reconstruction<>::output_type output(n_pixels * n_pixels, 0.0);
  output = r.emt();

  if (!cl.exist("output")) return 0;

  auto fn         = cl.get<cmdline::string>("output");
  auto fn_sep     = fn.find_last_of("\\/");
  auto fn_ext     = fn.find_last_of(".");
  auto fn_wo_ext  = fn.substr(0, fn_ext != std::string::npos
                             && (fn_sep == std::string::npos || fn_sep < fn_ext)
                               ? fn_ext : std::string::npos);
  std::ofstream out(fn);

  png_writer png(fn_wo_ext+".png");
  png.write_header<>(n_pixels, n_pixels);

  double output_max = 0.0;
  for(auto it = output.begin(); it != output.end(); ++it) {
    output_max = std::max(output_max, *it);
  }
  auto output_gain = static_cast<double>(std::numeric_limits<uint8_t>::max()) / output_max;

  size_t x = 0, y = 0;
  uint8_t row[n_pixels];
  for(auto it = output.begin(); it != output.end(); ++it) {

    if(y % n_pixels == 0 && y != 0) {
      x++; y = 0;
      png.write_row(row);
    }

    out << x << " " << y << " " << *it << std::endl;

    row[y] = std::numeric_limits<uint8_t>::max() - output_gain * *it;
    ++y;
  }
  png.write_row(row);

  return 0;

} catch(std::string &ex) {
  std::cerr << "error: " << ex << std::endl;
} catch(const char *ex) {
  std::cerr << "error: " << ex << std::endl;
}

}
