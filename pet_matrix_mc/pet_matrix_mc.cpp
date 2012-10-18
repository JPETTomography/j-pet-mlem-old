#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <vector>

#include <cmdline.h>

#include "event.h"
#include "detector.h"
#include "reconstruction.h"

int main(int argc, char *argv[]) {

  cmdline::parser cl;

  cl.add<int>("n-pixels",    'n', "number of pixels in one dimension", true);
  cl.add<int>("n-detectors", 'd', "number of ring detectors",          true);
  cl.add<int>("n-emissions", 'e', "emissions per pixel",               false, 1000000);

  cl.parse_check(argc, argv);

  std::cout << cl.get<int>("n-detectors") << " ";
  std::cout << cl.get<int>("n-pixels") << std::endl;

  int n_emissions = cl.get<int>("n-emissions");

  Ring2DDetector<int, double> detector(
    cl.get<int>("n-detectors"),
    cl.get<int>("n-pixels"));

  ProbabilityMatrix<int, double> probability_matrix(detector);

  const int max_lors = detector.n_detectors()*detector.n_detectors();

  double *p = new double[max_lors];

  for (int ipix = 0; ipix<probability_matrix.octant_size(); ++ipix) {
    auto pix = probability_matrix.octant(ipix);

    for (int i = 0; i<max_lors; ++i) p[i]=0.0;

    fill_pixel(pix, detector.n_detectors(), p, n_emissions);

    auto pixel_row = row_from_array(pix, p,detector.n_detectors(), max_lors);

    probability_matrix.push_back_row(pixel_row);
    std::cout << *pixel_row;
  }

  FILE *fout = fopen("prob.bin", "w");
  probability_matrix.fwrite(fout);
  fclose(fout);
}
