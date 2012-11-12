// PET Reconstruction
// Authors:
//   Adam Strzelecki <adam.strzlecki@uj.edu.pl>
//   Piotr Bialas <piotr.bialas@uj.edu.pl>
//
// Using Monte Carlo method and square detector scintilators

#include <cmdline.h>
#include "cmdline_types.h"

#if _OPENMP
#include <omp.h>
#endif

int main(int argc, char *argv[]) {

try {
  cmdline::parser cl;
  cl.footer("matrix_file measurements_file");

#if _OPENMP
  cl.add<size_t>     ("n-threads",   't', "number of OpenMP threads",          false);
#endif

  cl.parse_check(argc, argv);

#if _OPENMP
  if (cl.exist("n-threads")) {
    omp_set_num_threads(cl.get<size_t>("n-threads"));
  }
#endif

  // FIXME: IMPLEMENT ME!

  return 0;

} catch(std::string &ex) {
  std::cerr << "error: " << ex << std::endl;
} catch(const char *ex) {
  std::cerr << "error: " << ex << std::endl;
}

}
