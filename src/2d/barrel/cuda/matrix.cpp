#include <iostream>

#include "matrix.h"

#include "2d/barrel/detector_ring.h"
#include "2d/barrel/square_detector.h"
#include "2d/barrel/monte_carlo.h"
#include "2d/barrel/model.h"
#include "2d/geometry/point.h"

#include "util/progress.h"

namespace PET2D {
namespace Barrel {
namespace GPU {

OutputMatrix Matrix::run(cmdline::parser& cl) {

  DetectorRing detector_ring(cl.get<int>("n-detectors"),
                             cl.get<double>("radius"),
                             cl.get<double>("w-detector"),
                             cl.get<double>("h-detector"),
                             cl.get<double>("d-detector"));

  // GTX 770 - 8 SMX * 192 cores = 1536 cores -
  // each SMX can use 8 active blocks,

  auto n_blocks = cl.get<int>("cuda-blocks");
  auto n_threads_per_block = cl.get<int>("cuda-threads");
  auto n_threads = n_blocks * n_threads_per_block;

  auto& n_emissions = cl.get<int>("n-emissions");
  int n_emissions_per_thread = (n_emissions + n_threads - 1) / n_threads;

  // Number of emissions will be rounded to block size
  n_emissions = n_emissions_per_thread * n_blocks * n_threads_per_block;

  std::cout << "               blocks = " << n_blocks << std::endl
            << "    threads per block = " << n_threads_per_block << std::endl
            << " emissions per thread = " << n_emissions_per_thread << std::endl
            << "      total emissions = " << n_emissions << std::endl;

  auto tof_step = cl.get<double>("tof-step");
  int n_tof_positions = 1;
  double max_bias = 0;
  if (cl.exist("tof-step") && tof_step > 0) {
    max_bias = Model::max_bias();
    n_tof_positions = detector_ring.n_positions(tof_step, max_bias);
  }

  auto n_pixels = cl.get<int>("n-pixels");
  auto s_pixel = cl.get<double>("s-pixel");
  OutputMatrix output_matrix(
      n_pixels, cl.get<int>("n-detectors"), n_emissions, n_tof_positions);

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

  GPU::Matrix gpu_matrix(detector_ring,
                         n_threads_per_block,
                         n_blocks,
                         s_pixel,
                         n_tof_positions,
                         tof_step,
                         length_scale,
                         prng_seed.data());

  std::vector<int> pixel_hits(detector_ring.n_lors * n_tof_positions, 0);

  std::vector<LOR> lor_map;
  lor_map.resize(detector_ring.n_lors);
  for (LOR lor(0, 0); lor < LOR::end_for_detectors(detector_ring.n_detectors);
       ++lor) {
    lor_map[lor.index()] = lor;
  }

  auto end_pixel = Pixel::end_for_n_pixels_in_row(n_pixels / 2);
  auto n_pixels_total = end_pixel.index();
  util::progress progress(1, n_pixels_total);

  for (Pixel pixel(0, 0); pixel < end_pixel; ++pixel) {
    progress(pixel.index());

    gpu_matrix(pixel, n_emissions_per_thread, pixel_hits.data());

    for (int lor_index = 0; lor_index < lor_map.size(); ++lor_index) {
      auto lor = lor_map[lor_index];
      for (int position = 0; position < n_tof_positions; ++position) {
        auto hits = pixel_hits[n_tof_positions * lor_index + position];
        if (hits > 0) {
          output_matrix.emplace_back(lor, position, pixel, hits);
        }
      }
    }

    progress(pixel.index(), true);
  }

  return output_matrix;
}

}  // GPU
}  // Barrel
}  // PET2D
