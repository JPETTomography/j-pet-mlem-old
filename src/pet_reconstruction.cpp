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
#include "bstream.h"
#include "svg_ostream.h"
#include "cmdline_types.h"
#include "pet_rec.h"

#if _OPENMP
#include <omp.h>
#endif

std::vector< std::vector<hits_per_pixel> > Pet_reconstruction::system_matrix;
std::vector<int> Pet_reconstruction::list_of_lors;
std::vector<int> Pet_reconstruction::n;

int main(int argc, char *argv[]) {

try {
  cmdline::parser cl;

#if _OPENMP
  cl.add<size_t>     ("n-threads",'t', "number of OpenMP threads",          false);
#endif
  cl.add<std::string>("bin",'f', "system matrix binary file", false,"pet10k_full.bin");//"M_p128_d468_x100.bin");
  cl.add<std::string>("mean",'m', "mean file", false,"n_kuba.txt");
  cl.add("iterations",'n', "number of iterations", false,10);

  cl.parse_check(argc, argv);

#if _OPENMP
  if (cl.exist("n-threads")) {
    omp_set_num_threads(cl.get<size_t>("n-threads"));
  }
#endif

  Pet_reconstruction *pet_r = new Pet_reconstruction(cl.get<int>("iterations"),cl.get<std::string>("bin"),cl.get<std::string>("mean"));

  std::vector<float> output(pet_r->get_n_pixels() * pet_r->get_n_pixels(),0.0);

  output = pet_r->emt();

  ofstream data("out.txt");

  int x = 0;
  int y = 0;

  for(std::vector<float>::iterator it = output.begin(); it != output.end();++it)
  {

    if(y%128 == 0 && y != 0){x++; y = 0;}

    //printf("%d   %d   %f\n",x,y,*it);
    data << x << " " << y << " " << *it << endl;
    y++;

  }

  delete pet_r;

  // FIXME: IMPLEMENT ME!

  return 0;

} catch(std::string &ex) {
  std::cerr << "error: " << ex << std::endl;
} catch(const char *ex) {
  std::cerr << "error: " << ex << std::endl;
}

}
