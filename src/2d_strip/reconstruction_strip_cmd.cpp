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

#include "flags.h"
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

#endif
    cl.add("cpu", 'c', "run on cpu (via OPENMP)");
    cl.add<float>(
        "r-distance", 'r', "R distance between scientilators", false, 500.0f);
    cl.add<float>("s-length", 'l', "Scentilator_length", false, 1000.0f);
    cl.add<float>("p-size", 'p', "Pixel size", false, 5.0f);
    cl.add<int>("n-pixels", 'n', "Number of pixels", false, 200);
    cl.add<int>("iter", 'i', "number of iterations", false, 1);
    cl.add<float>("s-z", 's', "Sigma z error", false, 10.0f);
    cl.add<float>("s-dl", 'd', "Sigma dl error", false, 63.0f);
    cl.add<float>("gm", 'u', "Gamma error", false, 0.f);

    cl.parse_check(argc, argv);

#if _OPENMP
    if (cl.exist("n-threads")) {
      omp_set_num_threads(cl.get<int>("n-threads"));
    }
#endif

#if HAVE_CUDA
    if (cl.exist("gpu")) {

      float R_distance = cl.get<float>("r-distance");
      float Scentilator_length = cl.get<float>("s-length");
      float pixel_size = cl.get<float>("p-size");
      float n_pixels = Scentilator_length / pixel_size;
      float sigma = cl.get<float>("s-z");
      float dl = cl.get<float>("s-dl");

      Reconstruction<float> reconstruction(
          10, R_distance, Scentilator_length, n_pixels, pixel_size, sigma, dl);
      ibstream in("phantom.bin");

      in >> reconstruction;

      std::vector<event<float>>* event_data = reconstruction.get_data();

      std::cout << "VECTOR SIZE: " << event_data->size() << std::endl;

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
      cfg.inv_pow_sigma_dl = 1.0f / (dl * dl);
      cfg.inv_pow_sigma_z = 1.0f / (sigma * sigma);
      cfg.grid_size_y_ = n_pixels * pixel_size;
      cfg.grid_size_z_ = n_pixels * pixel_size;

      execute_kernel_reconstruction(
          cfg, event_data->data(), event_data->size());
    }
#endif

    if (cl.exist("cpu")) {
      float R_distance = cl.get<float>("r-distance");
      float Scentilator_length = cl.get<float>("s-length");
      float pixel_size = cl.get<float>("p-size");
      float n_pixels = Scentilator_length / pixel_size;
      float sigma = cl.get<float>("s-z");
      float dl = cl.get<float>("s-dl");

      int n_blocks = cl.get<int>("iter");
      Reconstruction<float> reconstruction(n_blocks,
                                           R_distance,
                                           Scentilator_length,
                                           n_pixels,
                                           pixel_size,
                                           sigma,
                                           dl);
      ibstream in("phantom.bin");

      in >> reconstruction;

      clock_t begin = clock();

      reconstruction(n_blocks);

      clock_t end = clock();

      std::cout << "Time:" << float(end - begin) / CLOCKS_PER_SEC / 4
                << std::endl;
    }
  } catch (std::string& ex) {
    std::cerr << "error: " << ex << std::endl;
  } catch (const char* ex) {
    std::cerr << "error: " << ex << std::endl;
  }

  return 0;
}
