#include <iostream>
#include "cmdline.h"
#include "util/cmdline_types.h"
#include "config.h"

// reference stuff from kernel.cu file
void run_kernel(char* str, int* val, int str_size, int val_size);

void phantom_kernel(int number_of_threads_per_block,
                    int blocks,
                    int n_emissions,
                    int pixels_in_row,
                    float radius,
                    float h_detector,
                    float w_detector,
                    float pixel_size);

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
    cl.add<float>("pixel-size", 's', "pixel size", false, 1.0f);
    cl.add<int>("block-num", 'B', "number of block", false, 64);
    cl.add<int>(
        "threads-per-block", 'P', "number of threads per block", false, 512);

    cl.parse_check(argc, argv);

    int pixels_in_row = cl.get<int>("n-pixels");
    int n_detectors = cl.get<int>("n-detectors");
    int n_emissions = cl.get<int>("n-emissions");
    float radius = cl.get<float>("radius");
    int s_pixel = cl.get<float>("s-pixel");
    float w_detector = cl.get<float>("w-detector");
    float h_detector = cl.get<float>("h-detector");
    // float tof_step = cl.get<float>("tof-step");

    int number_of_blocks = cl.get<int>("block-num");
    int number_of_threads_per_block = cl.get<int>("threads-per-block");

    char str[] = "Hello \0\0\0\0\0\0";
    int val[] = { 15, 10, 6, 0, -11, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    printf("%s\n", str);
    run_kernel(str, val, sizeof(str), sizeof(val));
    printf("%s\n", str);

    phantom_kernel(number_of_threads_per_block,
                   number_of_blocks,
                   n_emissions,
                   pixels_in_row,
                   radius,
                   h_detector,
                   w_detector,
                   s_pixel);
  }
  catch (std::string& ex) {

    std::cerr << "error: " << ex << std::endl;
  }
  catch (const char* ex) {

    std::cerr << "error: " << ex << std::endl;
  }

  return 0;
}
