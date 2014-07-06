#include <iostream>
#include <utility>
#include <iomanip>
#include <ctime>
#include <random>
#include <vector>
#include <exception>

#include "config.h"
#include "matrix.h"
#include "geometry.h"

#include "2d_xy/detector_ring.h"
#include "2d_xy/square_detector.h"
#include "2d_xy/monte_carlo.h"
#include "2d_xy/model.h"
#include "geometry/point.h"

double getwtime_cpu() {
#if !_MSC_VER
  struct timeval tv;
  static time_t sec = 0;
  gettimeofday(&tv, NULL);
  if (!sec)
    sec = tv.tv_sec;
  return (double)(tv.tv_sec - sec) + (double)tv.tv_usec / 1e6;
#else
  return 0;
#endif
}

bool run_gpu_matrix(int pixel_id,
                    int n_tof_positions,
                    int number_of_threads_per_block,
                    int number_of_blocks,
                    int n_emissions,
                    float radius,
                    float h_detector,
                    float w_detector,
                    float pixel_size,
                    gpu::LOR* lookup_table_lors,
                    Pixel<>* lookup_table_pixels,
                    unsigned int* cpu_prng_seed,
                    gpu::MatrixElement* cpu_matrix,
                    gpu::MatrixElement* gpu_output);

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

OutputMatrix run_gpu_matrix(cmdline::parser& cl) {

  auto pixels_in_row = cl.get<int>("n-pixels");
  auto n_detectors = cl.get<int>("n-detectors");
  auto n_emissions = cl.get<int>("n-emissions");
  auto radius = cl.get<double>("radius");
  auto s_pixel = cl.get<double>("s-pixel");
  auto w_detector = cl.get<double>("w-detector");
  auto h_detector = cl.get<double>("h-detector");

  int* emission_adr = &cl.get<int>("n-emissions");

  // GTX 770 - 8 SMX * 192 cores = 1536 cores -
  // each SMX can use 8 active blocks,

  auto number_of_blocks = cl.get<int>("cuda-blocks");
  auto number_of_threads_per_block = cl.get<int>("cuda-threads");

  int iteration_per_thread =
      floor(n_emissions / (number_of_blocks * number_of_threads_per_block));

  printf("Gpu grid config:\n");
  printf("Number of blocks:= %d\n", number_of_blocks);
  printf("Number of threads per block:= %d\n", number_of_threads_per_block);
  printf("Number of emissions per thread:= %d\n", iteration_per_thread);

  *emission_adr =
      iteration_per_thread * number_of_blocks * number_of_threads_per_block;

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

#if NO_TOF > 0
  int n_tof_positions = 1;
#else
  int n_tof_positions =
      ((int)(ceil(2.0f * (2.0f * (radius + h_detector)) / 0.01f)) + 1) / 2 * 2;
#endif

#if USE_GPU_KERNEL_PARAMETERS
  gpu_kernel_parameters kernel_parameters;

  kernel_parameters.x = 0;
  kernel_parameters.y = 0;
  kernel_parameters.iteration = iteration_per_thread;
  kernel_parameters.tof_n_positions = n_tof_positions;

  kernel_parameters.radius = radius;
  kernel_parameters.h_detector = h_detector;
  kernel_parameters.w_detector = w_detector;
  kernel_parameters.pixel_size = s_pixel;
#endif

  int triangle_pixel_size = (pixels_in_row / 2 * (pixels_in_row / 2 + 1) / 2);

  std::vector<gpu::LOR> lookup_table_lors;
  lookup_table_lors.resize(LORS);

  std::vector<Pixel<>> lookup_table_pixel;
  lookup_table_pixel.resize(triangle_pixel_size);

  gpu::MatrixElement matrix_element;

  for (int lor_i = 0; lor_i < LORS; ++lor_i) {

    matrix_element.hit[lor_i] = 0.f;
  }

#if NO_TOF > 0

  std::vector<gpu::MatrixElement> gpu_vector_output;
  gpu_vector_output.resize(1);

  std::vector<gpu::MatrixElement> cpu_matrix;
  cpu_matrix.resize(number_of_blocks);

#else

  std::vector<gpu::MatrixElement> gpu_vector_output;
  gpu_vector_output.resize(n_tof_positions);

  std::vector<gpu::MatrixElement> cpu_matrix;
  cpu_matrix.resize(n_tof_positions * number_of_blocks);

#endif

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

  OutputMatrix output_matrix(
      pixels_in_row, n_detectors, *emission_adr, n_tof_positions);

  double fulltime = double();

  for (int pixel_i = 0; pixel_i < triangle_pixel_size; ++pixel_i) {

    try {

      std::fill(
          gpu_vector_output.begin(), gpu_vector_output.end(), matrix_element);
      std::fill(cpu_matrix.begin(), cpu_matrix.end(), matrix_element);
    } catch (std::exception e) {
      std::cout << "Fill error:" << e.what() << std::endl;
    }

    std::default_random_engine gen;
    std::uniform_int_distribution<unsigned int> dis(1024, 1000000);
    gen.seed(12344 + 10 * pixel_i);

    for (int i = 0; i < 4 * number_of_blocks * number_of_threads_per_block;
         ++i) {

      cpu_prng_seed[i] = dis(gen);
    }

    double t0 = getwtime_cpu();

    run_gpu_matrix(pixel_i,
                   n_tof_positions,
                   number_of_threads_per_block,
                   number_of_blocks,
                   iteration_per_thread,
                   radius,
                   h_detector,
                   w_detector,
                   s_pixel,
                   lookup_table_lors.data(),
                   lookup_table_pixel.data(),
                   cpu_prng_seed,
                   cpu_matrix.data(),
                   gpu_vector_output.data());

    double t1 = getwtime_cpu();

    fulltime += (t1 - t0);

#if NO_TOF > 0

    for (int lor_i = 0; lor_i < LORS; ++lor_i) {

      if (gpu_vector_output[0].hit[lor_i] > 0) {

        OutputMatrix::LOR lor(lookup_table_lors[lor_i].lor_a,
                              lookup_table_lors[lor_i].lor_b);

        auto pixel = lookup_table_pixel[pixel_i];

        OutputMatrix::Element element;
        element.lor = lor;
        element.pixel = pixel;
        element.hits = gpu_vector_output[0].hit[lor_i];
        element.position = 0;

        output_matrix.push_back(element);
      }
    }
#else

    for (int tof_i = 0; tof_i < n_tof_positions; ++tof_i) {
      for (int lor_i = 0; lor_i < LORS; ++lor_i) {

        if (gpu_vector_output[tof_i].hit[lor_i] > 0) {

          OutputMatrix::LOR lor(lookup_table_lors[lor_i].lor_a,
                                lookup_table_lors[lor_i].lor_b);

          auto pixel = lookup_table_pixel[pixel_i];

          OutputMatrix::Element element;
          element.lor = lor;
          element.pixel = pixel;
          element.hits = gpu_vector_output[tof_i].hit[lor_i];
          element.position = tof_i;

          output_matrix.push_back(element);
        }
      }
    }

#endif
  }

  std::cout << fulltime << " " << fulltime / triangle_pixel_size << " "
            << fulltime / triangle_pixel_size /
                   (iteration_per_thread * number_of_blocks *
                    number_of_threads_per_block) << std::endl;

  free(cpu_prng_seed);

  return output_matrix;
}
