#include <iostream>
#include <utility>
#include <iomanip>
#include <ctime>
#include <random>
#include <vector>

#include "cmdline.h"
#include "util/cmdline_types.h"
#include "../util/png_writer.h"
#include "config.h"
#include "data_structures.h"

#include "2d_xy/detector_ring.h"
#include "2d_xy/square_detector.h"
#include "geometry/point.h"
#include "2d_xy/monte_carlo.h"
#include "2d_xy/sparse_matrix.h"
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
                    std::vector<Pixel<> >& lookup_table_pixel,
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
    std::vector<Pixel<> > lookup_table_pixel;
    lookup_table_pixel.resize(triangle_pixel_size);

    // Matrix_Element triangle_matrix_output[triangular_matrix_size];

    for (int j = pixels_in_row / 2 - 1; j >= 0; --j) {
      for (int i = 0; i <= j; ++i) {

        Pixel<> pixel(i, j);
        lookup_table_pixel[pixel.index()] = pixel;
      }
    }






    typedef SparseMatrix<Pixel<>, LOR<>> MatrixImpl;
    MatrixImpl matrix(pixels_in_row, n_detectors, n_tof_positions);


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


      for (auto p : gpu_vector_output) {

        for (int i = 0; i < LORS; ++i) {

          if (p.hit[i] > 0) {

            SparseElement<LOR<>,int,Pixel<>,int> temp;

            LOR<> lor(lookup_table_lors[i].lor_a,lookup_table_lors[i].lor_b);

            auto pixel = lookup_table_pixel[i];

            //printf("LOR(%d,%d) %f\n",lor.first,lor.second, p.hit[i]);

            temp.lor = lor;
            temp.pixel = pixel;
            temp.hits = p.hit[i];
            temp.position = 1;

            matrix.push_back(temp);
          }
        }

        obstream out("output.bin", std::ios::binary | std::ios::trunc);
        out << matrix;


#ifdef __TRASH__
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
#endif
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
