#include <iostream>
#include <utility>
#include <iomanip>
#include <ctime>
#include <random>
#include <vector>

#include "matrix_cuda.h"

#include "config.h"
#include "data_structures.h"

#include "2d_xy/detector_ring.h"
#include "2d_xy/square_detector.h"
#include "2d_xy/monte_carlo.h"
#include "2d_xy/model.h"
#include "geometry/point.h"

// reference stuff from kernel.cu file
// TODO: Should go into a separate header file.

void phantom_kernel(int number_of_threads_per_block,
                    int blocks,
                    int n_emissions,
                    int pixels_in_row,
                    float radius,
                    float h_detector,
                    float w_detector,
                    float pixel_size,
                    std::vector<Pixel<>>& lookup_table_pixel,
                    std::vector<Lor>& lookup_table_lors,
                    std::vector<Matrix_Element>& gpu_output);

void gpu_detector_geometry_kernel_test(float radius,
                                       float h_detector,
                                       float w_detector,
                                       float pixel_size,
                                       Detector_Ring& cpu_output);

void gpu_detector_hits_kernel_test(float crx,
                                   float cry,
                                   float cangle,
                                   float radius,
                                   float h_detector,
                                   float w_detector);

OutputMatrix run_gpu(cmdline::parser& cl) {

  auto pixels_in_row = cl.get<int>("n-pixels");
  auto n_detectors = cl.get<int>("n-detectors");
  auto n_emissions = cl.get<int>("n-emissions");
  auto radius = cl.get<double>("radius");
  auto s_pixel = cl.get<double>("s-pixel");
  auto w_detector = cl.get<double>("w-detector");
  auto h_detector = cl.get<double>("h-detector");

  auto number_of_blocks = cl.get<int>("n-blocks") ?: 512;
  auto number_of_threads_per_block = cl.get<int>("n-threads");

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

  // Lor lookup_table_lors[LORS];
  std::vector<Lor> lookup_table_lors;
  lookup_table_lors.resize(LORS);

  for (int i = 0; i < NUMBER_OF_DETECTORS; ++i) {
    for (int j = 0; j < NUMBER_OF_DETECTORS; ++j) {

      Lor temp;
      temp.lor_a = i;
      temp.lor_b = j;
      lookup_table_lors[(i * (i + 1) / 2) + j] = temp;
    }
  }

  // Pixel<> lookup_table_pixel[triangular_matrix_size];
  std::vector<Pixel<>> lookup_table_pixel;
  lookup_table_pixel.resize(triangle_pixel_size);

  // Matrix_Element triangle_matrix_output[triangular_matrix_size];

  for (int j = pixels_in_row / 2 - 1; j >= 0; --j) {
    for (int i = 0; i <= j; ++i) {

      Pixel<> pixel(i, j);
      lookup_table_pixel[pixel.index()] = pixel;
    }
  }

  std::vector<Matrix_Element> gpu_vector_output;

  gpu_vector_output.resize(triangle_pixel_size);

  std::cout << "VECTOR size: " << gpu_vector_output.size() << std::endl;

  phantom_kernel(number_of_threads_per_block,
                 number_of_blocks,
                 n_emissions,
                 pixels_in_row,
                 radius,
                 h_detector,
                 w_detector,
                 s_pixel,
                 lookup_table_pixel,
                 lookup_table_lors,
                 gpu_vector_output);

  OutputMatrix output_matrix(pixels_in_row, n_detectors, n_tof_positions);

  for (auto p : gpu_vector_output) {

    for (int i = 0; i < LORS; ++i) {

      if (p.hit[i] > 0) {

        OutputMatrix::LOR lor(lookup_table_lors[i].lor_a, lookup_table_lors[i].lor_b);

        auto pixel = lookup_table_pixel[i];

        // printf("LOR(%d,%d) %f\n",lor.first,lor.second, p.hit[i]);

        OutputMatrix::Element element;
        element.lor = lor;
        element.pixel = pixel;
        element.hits = p.hit[i];
        element.position = 0;

        output_matrix.push_back(element);
      }
    }
  }

  return output_matrix;
}
