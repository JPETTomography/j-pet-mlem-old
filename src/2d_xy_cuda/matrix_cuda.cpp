#include <iostream>
#include <utility>
#include <iomanip>
#include <ctime>
#include <random>
#include <vector>

#include "cmdline.h"
#include "util/cmdline_types.h"
#include "config.h"
#include "data_structures.h"

#include "2d_xy/detector_ring.h"
#include "2d_xy/square_detector.h"
#include "geometry/point.h"
#include "2d_xy/monte_carlo.h"
#include "2d_xy/matrix_pixel_major.h"
#include "2d_xy/model.h"
#include "2d_xy/lor.h"
#include "geometry/pixel.h"

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

int main(int argc, char* argv[]) {

  try {

    cmdline::parser cl;
    cl.footer("matrix_file ...");

    cl.add<cmdline::string>(
        "config", 'c', "load config file", false, cmdline::string(), false);
    cl.add<int>(
        "n-pixels", 'n', "number of pixels in one dimension", false, 64);
    cl.add<int>("n-detectors", 'd', "number of ring detectors", false, 64);
    cl.add<int>("n-emissions", 'e', "emissions per pixel", false, 1000);
    cl.add<float>("radius", 'r', "inner detector ring radius", false, 100);
    cl.add<float>("s-pixel", 'p', "pixel size", false, 1.0f);
    cl.add<float>("tof-step", 'T', "TOF quantisation step", false);
    cl.add<float>("w-detector", 'w', "detector width", false, 1.0f);
    cl.add<float>("h-detector", 'h', "detector height", false, 2.0f);
    cl.add<int>("block-num", 'B', "number of block", false, 64);
    cl.add<int>(
        "threads-per-block", 'P', "number of threads per block", false, 512);

    cl.parse_check(argc, argv);

    int pixels_in_row = cl.get<int>("n-pixels");
    int n_detectors = cl.get<int>("n-detectors");
    int n_emissions = cl.get<int>("n-emissions");
    float radius = cl.get<float>("radius");
    float s_pixel = cl.get<float>("s-pixel");
    float w_detector = cl.get<float>("w-detector");
    float h_detector = cl.get<float>("h-detector");
    float tof_step = 0;

    int number_of_blocks = cl.get<int>("block-num");
    int number_of_threads_per_block = cl.get<int>("threads-per-block");

    std::cerr << radius << " " << s_pixel << " " << w_detector << " "
              << h_detector << std::endl;

    int n_tof_positions = 1;

    typedef MatrixPixelMajor<Pixel<>, LOR<>> MatrixImpl;
    MatrixImpl matrix(pixels_in_row, n_detectors, n_tof_positions);

    std::vector<Matrix_Element> gpu_vector_output;

    gpu_vector_output.resize(matrix.total_n_pixels_in_triangle());

    std::cout << "VECTOR size: " << gpu_vector_output.size() << std::endl;

    phantom_kernel(number_of_threads_per_block,
                   number_of_blocks,
                   n_emissions,
                   pixels_in_row,
                   radius,
                   h_detector,
                   w_detector,
                   s_pixel,
                   gpu_vector_output);

    DetectorRing<> dr(n_detectors, radius, w_detector, h_detector);

    MonteCarlo<DetectorRing<>, MatrixImpl> monte_carlo(
        dr, matrix, s_pixel, tof_step);

    std::random_device rd;
    tausworthe gen(rd());
    gen.seed(2345255);

    long cpu_emissions =
        (long)n_emissions * number_of_blocks * number_of_threads_per_block;

    std::cerr << "CPU " << cpu_emissions << std::endl;

    clock_t begin = clock();

    monte_carlo.emit_pixel(gen, AlwaysAccept<>(), cpu_emissions);

    matrix.get_pixel_data(cpu_emissions, 0);

    clock_t end = clock();

    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    std::cout << "CPU TIME: " << elapsed_secs << std::endl;

    std::cout << "KERNEL OUTPUT" << std::endl;

    for (auto p : gpu_vector_output) {

      for (int i = 0; i < LORS; ++i) {

        if (p.hit[i] > 0) {

          printf("%f\n", p.hit[i]);
        }
      }
    }
  }

  catch (std::string& ex) {

    std::cerr << "error: " << ex << std::endl;
  }
  catch (const char* ex) {

    std::cerr << "error: " << ex << std::endl;
  }

  return 0;
}
