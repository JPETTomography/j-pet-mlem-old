#include <iostream>

#include "matrix.h"

#include "3d/hybrid/options.h"

#include "2d/barrel/ring_scanner.h"
#include "2d/barrel/square_detector.h"
#include "2d/barrel/monte_carlo.h"
#include "common/model.h"
#include "2d/geometry/point.h"
#include "2d/barrel/scanner_builder.h"

#include "util/progress.h"

namespace PET3D {
namespace Hybrid {
namespace GPU {

/// \cond PRIVATE

OutputMatrix Matrix::run(cmdline::parser& cl) {

  Scanner scanner(
      PET2D::Barrel::ScannerBuilder<Scanner2D>::build_multiple_rings(
          PET3D_LONGITUDINAL_SCANNER_CL(cl, Scanner::F)),
      F(cl.get<double>("length")));

  // GTX 770 - 8 SMX * 192 cores = 1536 cores -
  // each SMX can use 8 active blocks,

  auto n_blocks = cl.get<int>("cuda-blocks");
  auto n_threads_per_block = cl.get<int>("cuda-threads");
  auto n_threads = n_blocks * n_threads_per_block;

  auto& n_emissions = cl.get<int>("n-emissions");
  int n_thread_emissions = (n_emissions + n_threads - 1) / n_threads;

  // Number of emissions will be rounded to block size
  n_emissions = n_thread_emissions * n_blocks * n_threads_per_block;

  if (cl.exist("verbose")) {
    std::cout << "GPU parameters:" << std::endl;
    std::cout << "           blocks = " << n_blocks << std::endl
              << "threads per block = " << n_threads_per_block << std::endl
              << " thread emissions = " << n_thread_emissions << std::endl
              << "  total emissions = " << n_emissions << std::endl;
  }

  auto z_position = cl.get<double>("z-position");
  auto n_pixels = cl.get<int>("n-pixels");
  auto s_pixel = cl.get<double>("s-pixel");
  OutputMatrix output_matrix(
      n_pixels, cl.get<int>("n-detectors"), n_emissions, 1);

  // convert obsolete acceptance option to length scale
  auto& length_scale = cl.get<double>("base-length");
  if (cl.exist("acceptance") && !cl.exist("base-length")) {
    length_scale = 1.0 / cl.get<double>("acceptance");
  }

  std::vector<unsigned int> prng_seed(n_blocks * n_threads_per_block * 4);

  util::random::tausworthe gen;
  gen.seed(345555);
  for (int i = 0; i < 4 * n_blocks * n_threads_per_block; ++i) {
    prng_seed[i] = gen();  // 53445 + i
  }

  GPU::Matrix gpu_matrix(z_position,
                         scanner,
                         n_threads_per_block,
                         n_blocks,
                         s_pixel,
                         length_scale,
                         prng_seed.data());

  const auto n_lors = LOR::end_for_detectors(scanner.barrel.size()).index();
  std::vector<int> pixel_hits(n_lors, 0);

  std::vector<LOR> lor_map;
  lor_map.resize(n_lors);
  for (LOR lor(0, 0); lor < LOR::end_for_detectors(scanner.barrel.size());
       ++lor) {
    lor_map[lor.index()] = lor;
  }

  auto end_pixel = Pixel::end_for_n_pixels_in_row(n_pixels / 2);
  auto n_pixels_total = end_pixel.index();
  util::progress progress(1, n_pixels_total);

  for (Pixel pixel(0, 0); pixel < end_pixel; ++pixel) {
    progress(pixel.index());

    gpu_matrix(pixel, n_thread_emissions, pixel_hits.data());

    for (size_t lor_index = 0; lor_index < lor_map.size(); ++lor_index) {
      auto lor = lor_map[lor_index];
      auto hits = pixel_hits[lor_index];
      if (hits > 0) {
        output_matrix.emplace_back(lor, 0, pixel, hits);
      }
    }

    progress(pixel.index(), true);
  }

  return output_matrix;
}

/// \endcond

}  // GPU
}  // Barrel
}  // PET2D
