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
#include "reconstruction.h"
#include "png_writer.h"

#if _OPENMP
#include <omp.h>
#endif

template<typename Iterator> 
void write_text_from_vector(std::ostream &out,Iterator start, Iterator stop, int line_length) {
  int index=0;
  for(;start!=stop;++start,++index) {
    out<<*start<<" ";
    if(index%line_length==0)
      out<<"\n";
  }
}


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

  reconstruction<> reconstructor(
    cl.get<size_t>("iterations"),
    cl.get<cmdline::string>("system"),
    cl.get<cmdline::string>("mean"));

  auto n_pixels = reconstructor.get_n_pixels();
  auto total_n_pixels=n_pixels*n_pixels;
  reconstruction<>::output_type rho(total_n_pixels, 0.0);
  reconstruction<>::output_type rho_detected(total_n_pixels, 0.0);
  reconstructor.emt(cl.get<size_t>("iterations"));
  rho = reconstructor.rho();
  rho_detected = reconstructor.rho_detected();

  if (!cl.exist("output")) return 0;

  auto fn         = cl.get<cmdline::string>("output");
  auto it_fn_sep     = fn.find_last_of("\\/");
  auto it_fn_ext     = fn.find_last_of(".");
  auto fn_wo_ext  = fn.substr(0, it_fn_ext != std::string::npos
                             && (it_fn_sep == std::string::npos || it_fn_sep < it_fn_ext)
                               ? it_fn_ext : std::string::npos);
  auto fn_ext     = fn.substr(it_fn_ext != std::string::npos?it_fn_ext:fn.size(),fn.size());

  std::ofstream out(fn);
  write_text_from_vector(out,rho.begin(),rho.end(),n_pixels);
  out.close();

  std::ofstream out_detected(fn_wo_ext+"_detected"+fn_ext);
  write_text_from_vector(out_detected,rho_detected.begin(),rho_detected.end(),n_pixels);
  out_detected.close();


  // output reconstruction PNG
  png_writer png(fn_wo_ext+".png");
  png.write_header<>(n_pixels, n_pixels);


  double output_max = 0.0;
  for(auto it = rho.begin(); it != rho.end(); ++it) {
    output_max = std::max(output_max, *it);
  }

  auto output_gain = static_cast<double>(std::numeric_limits<uint8_t>::max()) / output_max;

  for (int y = n_pixels-1; y >= 0; --y) {
    uint8_t row[n_pixels];
    for (auto x = 0; x < n_pixels; ++x) {
      row[x] = std::numeric_limits<uint8_t>::max() - output_gain * rho[LOCATION(x,y,n_pixels)];
    }
    png.write_row(row);
  }

  return 0;

} catch(std::string &ex) {
  std::cerr << "error: " << ex << std::endl;
} catch(const char *ex) {
  std::cerr << "error: " << ex << std::endl;
}

}
