#include <iostream>
#include <utility>
#include <iomanip>
#include <ctime>
#include <random>

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

void phantom_kernel(int number_of_threads_per_block,
                    int blocks,
                    int n_emissions,
                    int pixels_in_row,
                    float radius,
                    float h_detector,
                    float w_detector,
                    float pixel_size);

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
    cl.add<int>("n-emissions", 'e', "emissions per pixel", false, 10000);
    cl.add<float>("radius", 'r', "inner detector ring radius", false, 100);
    cl.add<float>("s-pixel", 'p', "pixel size", false, 1.0f);
    cl.add<float>("tof-step", 'T', "TOF quantisation step", false);
    cl.add<float>("w-detector", 'w', "detector width", false, 1.0f);
    cl.add<float>("h-detector", 'h', "detector height", false, 1.0f);
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

    std::cerr<<radius<<" "<<s_pixel<<" "<<w_detector<<" "<<h_detector<<std::endl;
    //----------SYSTEM MATRIX GENERATION----------//

    phantom_kernel(number_of_threads_per_block,
                   number_of_blocks,
                   n_emissions,
                   pixels_in_row,
                   radius,
                   h_detector,
                   w_detector,
                   s_pixel);

    DetectorRing<> dr(n_detectors, radius, w_detector, h_detector);

#ifdef TEST

    Detector_Ring cpu_output;

    typedef std::pair<int, float> error;

    std::vector<error> error_list;

    float epsilon_error = 0.0001f;

    gpu_detector_geometry_kernel_test(
        radius, h_detector, w_detector, s_pixel, cpu_output);

    for (int detector_id = 0; detector_id < NUMBER_OF_DETECTORS;
         ++detector_id) {

      auto detector_points = dr[detector_id];

#if VERBOSE
      std::cout << "DETECTOR: " << detector_id << std::endl;
#endif

      for (int point_id = 0; point_id < 4; ++point_id) {

        auto point = detector_points[point_id];
        float diff = std::fabs(
            point.x - cpu_output.detector_list[detector_id].points[point_id].x);
        if (diff > epsilon_error) {
          error_list.push_back(std::make_pair(detector_id, diff));
#if VERBOSE
          std::cout << "Diff x : " << diff << std::endl;
#endif
        }

        diff = std::fabs(
            point.y - cpu_output.detector_list[detector_id].points[point_id].y);
        if (diff > epsilon_error) {
          error_list.push_back(std::make_pair(detector_id, diff));
#if VERBOSE
          std::cout << "Diff y : " << diff << std::endl;
#endif
        }

#if VERBOSE

        std::cout << std::setprecision(10) << "Cpu representation: " << point.x
                  << " " << point.y << std::endl;
        std::cout << std::setprecision(10) << "Gpu representation: "
                  << cpu_output.detector_list[detector_id].points[point_id].x
                  << " "
                  << cpu_output.detector_list[detector_id].points[point_id].y
                  << std::endl;
#endif
      }
    }

    if (!error_list.size()) {

      std::cout << "Number of errors in cpu|gpu detectors geometry comparison: "
                << error_list.size() << std::endl;
    }

    //----------MATRIX OUTPUT----------//

    std::cout << "Matrix output Test:" << std::endl;

#endif

    int n_tof_positions = 1;

    typedef MatrixPixelMajor<Pixel<>, LOR<>> MatrixImpl;
    MatrixImpl matrix(pixels_in_row, n_detectors, n_tof_positions);
    MonteCarlo<DetectorRing<>, MatrixImpl> monte_carlo(
        dr, matrix, s_pixel, tof_step);

    std::random_device rd;
    tausworthe gen(rd());
    gen.seed(2345255);

    //    HIT DATA:0.149976 0.59646 0.979102

    int cpu_emmisions = 10000 * 512 * 64;

    //    monte_carlo.emit_pixel(
    //        gen, AlwaysAccept<>(), cpu_emmisions, 0.149976, 0.59646,
    // 0.979102);

    clock_t begin = clock();

    monte_carlo.emit_pixel(gen, AlwaysAccept<>(), cpu_emmisions);

    matrix.get_pixel_data(cpu_emmisions, 1);

    clock_t end = clock();

    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    std::cout << "CPU TIME: " << elapsed_secs << std::endl;

#ifdef TEST

    //-------------GPU_DETECTOR_HITS_TEST------------------//

    gpu_detector_hits_kernel_test(
        0.149976, 0.59646, 0.979102, radius, h_detector, w_detector);

#endif
  }

  catch (std::string& ex) {

    std::cerr << "error: " << ex << std::endl;
  }
  catch (const char* ex) {

    std::cerr << "error: " << ex << std::endl;
  }

  return 0;
}
