#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <vector>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
namespace po = boost::program_options;

#include "event.h"
#include "detector.h"
#include "reconstruction.h"

int main(int argc, char *argv[]) {

  int n_emissions;
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",                          "shows this help")
    ("n-pixels,n",    po::value<int>(), "set number of pixels in one dimension")
    ("n-detectors,d", po::value<int>(), "set number of ring detector")
    ("n-emissions,e", po::value<int>(&n_emissions)->default_value(1000000), "emissions per pixel")
  ;

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
  } catch(const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 1;
  }

  if (!vm.count("n-detectors") || !vm.count("n-pixels")) {
    std::cerr << "missing number of detectors or number of pixels" << std::endl;
    return 1;
  }

  std::cout << vm["n-detectors"].as<int>() << " ";
  std::cout << vm["n-pixels"].as<int>() << std::endl;

  Ring2DDetector<int, double> detector(
    vm["n-detectors"].as<int>(),
    vm["n-pixels"].as<int>());

  ProbabilityMatrix<int, double> probability_matrix(detector);

  const int max_lors = detector.n_detectors()*detector.n_detectors();

  double *p = new double[max_lors];

  for (int ipix = 0; ipix<probability_matrix.octant_size(); ++ipix) {
    auto pix = probability_matrix.octant(ipix);

    for (int i = 0; i<max_lors; ++i) p[i]=0.0;

    fill_pixel(pix, detector.n_detectors(), p,n_emissions);

    auto pixel_row = row_from_array(pix, p,detector.n_detectors(), max_lors);

    probability_matrix.push_back_row(pixel_row);
    std::cout << *pixel_row;
  }

  FILE *fout = fopen("prob.bin", "w");
  probability_matrix.fwrite(fout);
  fclose(fout);
}
