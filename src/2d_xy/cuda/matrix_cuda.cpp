#include <iostream>
#include <utility>
#include <iomanip>
#include <ctime>
#include <random>
#include <vector>

#include "config.h"
#include "matrix.h"
#include "geometry.h"

#include "2d_xy/detector_ring.h"
#include "2d_xy/square_detector.h"
#include "2d_xy/monte_carlo.h"
#include "2d_xy/model.h"
#include "geometry/point.h"

#if FIXME
// reference stuff from kernel.cu file
// TODO: Should go into a separate header file.

void run_phantom_kernel(int number_of_threads_per_block,
                        int blocks,
                        int n_emissions,
                        int pixels_in_row,
                        float radius,
                        float h_detector,
                        float w_detector,
                        float pixel_size,
                        Pixel<>* lookup_table_pixel,
                        Lor* lookup_table_lors,
                        std::vector<MatrixElement>& gpu_output);
#endif

void run_monte_carlo_kernel(int number_of_threads_per_block,
                            int number_of_blocks,
                            int n_emissions,
                            int n_detectors,
                            int pixels_in_row,
                            int triangle_pixel_size,
                            float radius,
                            float h_detector,
                            float w_detector,
                            float pixel_size,
                            gpu::LOR* lookup_table_lors,
                            Pixel<>* lookup_table_pixels,
                            unsigned int* cpu_prng_seed,
                            gpu::MatrixElement* cpu_matrix,
                            gpu::MatrixElement* gpu_output);

void run_detector_geometry_test_kernel(float radius,
                                       float h_detector,
                                       float w_detector,
                                       float pixel_size,
                                       gpu::DetectorRing& cpu_output);

void run_detector_hits_test_kernel(float crx,
                                   float cry,
                                   float cangle,
                                   float radius,
                                   float h_detector,
                                   float w_detector);

void fill_gpu_data(gpu::LOR* lookup_table_lors,
                   Pixel<>* lookup_table_pixel,
                   unsigned int* cpu_prng_seed,
                   int& number_of_blocks,
                   int& number_of_threads_per_block,
                   int& pixels_in_row) {

  for (int i = 0; i < NUMBER_OF_DETECTORS; ++i) {
    for (int j = 0; j < NUMBER_OF_DETECTORS; ++j) {

      gpu::LOR temp;
      temp.lor_a = i;
      temp.lor_b = j;
      lookup_table_lors[(i * (i + 1) / 2) + j] = temp;
    }
  }

  for (int j = pixels_in_row / 2 - 1; j >= 0; --j) {
    for (int i = 0; i <= j; ++i) {

      Pixel<> pixel(i, j);
      lookup_table_pixel[pixel.index()] = pixel;
    }
  }

  std::default_random_engine gen;
  std::uniform_int_distribution<unsigned int> dis(1024, 1000000);
  gen.seed(345555);

  for (int i = 0; i < 4 * number_of_blocks * number_of_threads_per_block; ++i) {

    // cpu_prng_seed[i] = 53445 + i; //dis(gen);
    cpu_prng_seed[i] = dis(gen);
  }
}

OutputMatrix run_gpu(cmdline::parser& cl) {

  auto pixels_in_row = cl.get<int>("n-pixels");
  auto n_detectors = cl.get<int>("n-detectors");
  auto n_emissions = cl.get<int>("n-emissions");
  auto radius = cl.get<double>("radius");
  auto s_pixel = cl.get<double>("s-pixel");
  auto w_detector = cl.get<double>("w-detector");
  auto h_detector = cl.get<double>("h-detector");



  // GTX 770 - 8 SMX * 192 cores = 1536 cores -
  // each SMX can use 8 active blocks,



  auto number_of_blocks = cl.get<int>("n-blocks") ?: 96;
  auto number_of_threads_per_block = cl.get<int>("n-threads") ?: 512;

  auto iteration_per_thread = floor(n_emissions/ (number_of_blocks * number_of_threads_per_block));

  // automatic pixel size
  if (!cl.exist("radius")) {
    if (!cl.exist("s-pixel")) {
      radius = M_SQRT2;  // exact result
    } else {
      radius = s_pixel * pixels_in_row / M_SQRT2;
    }
    std::cerr << "--radius=" << radius << std::endl;
  }

  // automatic radius
  if (!cl.exist("s-pixel")) {
    if (!cl.exist("radius")) {
      s_pixel = 2. / pixels_in_row;  // exact result
    } else {
      s_pixel = M_SQRT2 * radius / pixels_in_row;
    }
    std::cerr << "--s-pixel=" << s_pixel << std::endl;
  }

  // automatic detector size
  if (!cl.exist("w-detector")) {
    w_detector = 2 * M_PI * .9 * radius / n_detectors;
    std::cerr << "--w-detector=" << w_detector << std::endl;
  }
  if (!cl.exist("h-detector")) {
    h_detector = w_detector;
    std::cerr << "--h-detector=" << h_detector << std::endl;
  }

  std::cerr << radius << " " << s_pixel << " " << w_detector << " "
            << h_detector << std::endl;

  int n_tof_positions = 1;
  int triangle_pixel_size = (pixels_in_row / 2 * (pixels_in_row / 2 + 1) / 2);

  std::vector<gpu::LOR> lookup_table_lors;
  lookup_table_lors.resize(LORS);

  std::vector<Pixel<>> lookup_table_pixel;
  lookup_table_pixel.resize(triangle_pixel_size);

  std::vector<gpu::MatrixElement> gpu_vector_output;
  gpu_vector_output.resize(triangle_pixel_size);

  unsigned int* cpu_prng_seed;

  cpu_prng_seed =
      (unsigned int*)malloc(number_of_blocks * number_of_threads_per_block * 4 *
                            sizeof(unsigned int));

  fill_gpu_data(lookup_table_lors.data(),
                lookup_table_pixel.data(),
                cpu_prng_seed,
                number_of_blocks,
                number_of_threads_per_block,
                pixels_in_row);

  for (int i = 0; i < triangle_pixel_size; ++i) {

    for (int lor = 0; lor < LORS; ++lor) {

      gpu_vector_output[i].hit[lor] = 0;
    }
  }

  std::vector<gpu::MatrixElement> cpu_matrix;
  cpu_matrix.resize(number_of_blocks);

  run_monte_carlo_kernel(number_of_threads_per_block,
                         number_of_blocks,
                         iteration_per_thread,
                         n_detectors,
                         pixels_in_row,
                         triangle_pixel_size,
                         radius,
                         h_detector,
                         w_detector,
                         s_pixel,
                         lookup_table_lors.data(),
                         lookup_table_pixel.data(),
                         cpu_prng_seed,
                         cpu_matrix.data(),
                         gpu_vector_output.data());

  OutputMatrix output_matrix(pixels_in_row, n_detectors, n_tof_positions);

  for (int id = 0; id < triangle_pixel_size; ++id) {

    for (int i = 0; i < LORS; ++i) {

      if (gpu_vector_output[id].hit[i] > 0) {

        OutputMatrix::LOR lor(lookup_table_lors[i].lor_a,
                              lookup_table_lors[i].lor_b);

        auto pixel = lookup_table_pixel[id];
#ifdef PRINT
        if (id == 1) {
          printf("ID: %d LOR(%d,%d) %f\n",
                 id,
                 lor.first,
                 lor.second,
                 gpu_vector_output[id].hit[i]);
        }
#endif
        OutputMatrix::Element element;
        element.lor = lor;
        element.pixel = pixel;
        element.hits = gpu_vector_output[id].hit[i];
        element.position = 0;

        output_matrix.push_back(element);
      }
    }
  }

  free(cpu_prng_seed);

  return output_matrix;
}
