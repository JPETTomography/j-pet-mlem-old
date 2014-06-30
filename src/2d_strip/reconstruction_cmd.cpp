#include <iostream>
#include <vector>
#include <ctime>

#if SSE_FLUSH
#include <xmmintrin.h>
#endif

#include "cmdline.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"
#include "util/cmdline_types.h"
#include "util/png_writer.h"

#include "event.h"
#include "reconstruction.h"
#include "config.h"

#if HAVE_CUDA
#include "gpu_kernel_wrapper.h"
#endif

using namespace std;

int main(int argc, char* argv[]) {

#if SSE_FLUSH
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif

  try {
    cmdline::parser cl;

#if _OPENMP
    cl.add<int>("n-threads", 't', "number of OpenMP threads", false, 4);
#endif
#if HAVE_CUDA
    cl.add("gpu", 'g', "run on GPU (via CUDA)");
    cl.add<int>("cuda-blocks", 'b', "number of CUDA blocks", false, 1);
    cl.add<int>(
        "cuda-threads", 'w', "number of CUDA threads per block", false, 512);
    cl.add<int>("warp-offset", 0, "warp offset for test only", false, 1);

#endif
    cl.add<double>(
        "r-distance", 'r', "R distance between scientilators", false, 500);
    cl.add<double>("s-length", 'l', "Scentilator_length", false, 1000);
    cl.add<double>("p-size", 'p', "Pixel size", false, 5);
    cl.add<int>("n-pixels", 'n', "Number of pixels", false, 200);
    cl.add<int>("iter", 'i', "number of iterations", false, 1);
    cl.add<double>("s-z", 's', "Sigma z error", false, 10);
    cl.add<double>("s-dl", 'd', "Sigma dl error", false, 63);
    cl.add<double>("gm", 'u', "Gamma error", false, 0);
    cl.add<cmdline::string>(
        "output", 'o', "output files (png)", false, "cpu_rec_iteration");
    cl.parse_check(argc, argv);

#if _OPENMP
    if (cl.exist("n-threads")) {
      omp_set_num_threads(cl.get<int>("n-threads"));
    }
#endif

    double R_distance = cl.get<double>("r-distance");
    double Scentilator_length = cl.get<double>("s-length");
    double pixel_size = cl.get<double>("p-size");
    double n_pixels = Scentilator_length / pixel_size;
    double sigma = cl.get<double>("s-z");
    double dl = cl.get<double>("s-dl");

    Reconstruction<double> reconstruction(
        R_distance, Scentilator_length, n_pixels, pixel_size, sigma, dl);

    for (auto& fn : cl.rest()) {
      ibstream in(fn);
      reconstruction << in;
    }

#if HAVE_CUDA
    if (cl.exist("gpu")) {
      gpu_config::GPU_parameters cfg;
      cfg.R_distance = R_distance;
      cfg.Scentilator_length = Scentilator_length;
      cfg.pixel_size = pixel_size;
      cfg.n_pixels = n_pixels;
      cfg.sigma = sigma;
      cfg.dl = dl;
      cfg.number_of_blocks = cl.get<int>("cuda-blocks");
      cfg.number_of_threads_per_block = cl.get<int>("cuda-threads");
      cfg.number_of_events = 1;
      cfg.inv_pow_sigma_dl = 1 / (dl * dl);
      cfg.inv_pow_sigma_z = 1 / (sigma * sigma);
      cfg.grid_size_y_ = n_pixels * pixel_size;
      cfg.grid_size_z_ = n_pixels * pixel_size;

      execute_kernel_reconstruction(cfg,
                                    reconstruction.get_event_list(),
                                    cl.get<int>("warp-offset"),
                                    cl.get<int>("iter"));
    } else
#endif
    {
      clock_t begin = clock();
      reconstruction(cl.get<int>("iter"), 1, cl.get<cmdline::string>("output"));
      clock_t end = clock();

      std::cout << "Time: " << double(end - begin) / CLOCKS_PER_SEC / 4
                << std::endl;
    }
  } catch (std::string& ex) {
    std::cerr << "error: " << ex << std::endl;
  } catch (const char* ex) {
    std::cerr << "error: " << ex << std::endl;
  }

  return 0;
}
