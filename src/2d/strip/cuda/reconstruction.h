#pragma once

#include <vector>
#include <sstream>
#include <iomanip>

#include "util/png_writer.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"
#include "util/progress.h"

#include "../response.h"
#include "../scanner.h"

namespace PET2D {
namespace Strip {
/// CUDA optimized subimplementation
namespace GPU {

/// CUDA entry-point function
template <typename F>
void run_reconstruction(Scanner<F, short>& scanner,
                        Response<F>* responses,
                        int n_responses,
                        int n_iteration_blocks,
                        int n_iterations_in_block,
                        void (*output_callback)(Scanner<F, short>& scanner,
                                                int iteration,
                                                F* image,
                                                void* context),
                        void (*progress_callback)(int iteration,
                                                  void* context,
                                                  bool finished),
                        void* context,
                        int device,
                        int n_blocks,
                        int n_threads_per_block,
                        bool verbose);

struct Context {
  Context(util::progress& progress,
          std::string& output_file_name,
          bool text_output,
          int n_file_digits)
      : progress(progress),
        output_file_name(output_file_name),
        text_output(text_output),
        n_file_digits(n_file_digits) {}
  util::progress& progress;
  std::string& output_file_name;
  bool text_output;
  int n_file_digits;
};

void output(Scanner<float, short>& scanner,
            int iteration,
            float* output,
            void* ptr) {
  Context* context = static_cast<Context*>(ptr);

  if (!context->output_file_name.length()) {
    return;
  }

  std::stringstream base_name;
  base_name << context->output_file_name << "_";  // phantom_
  if (iteration >= 0) {
    base_name << std::setw(context->n_file_digits)  //
              << std::setfill('0')                  //
              << iteration << std::setw(0);         // 001
  } else {
    base_name << "sensitivity";
  }

  if (context->text_output) {
    std::ofstream txt(base_name.str() + ".txt");
    txt.precision(12);
    txt << std::fixed;
    for (int y = 0; y < scanner.n_y_pixels; ++y) {
      for (auto x = 0; x < scanner.n_z_pixels; ++x) {
        auto value = output[y * scanner.n_z_pixels + x];
        if (value >= 0.000000000001f) {
          txt << x << ' ' << y << ' ' << value << std::endl;
        }
      }
    }
  } else {
    util::obstream bin(base_name.str() + ".bin");
    bin.write(output, scanner.total_n_pixels);
  }

  util::png_writer png(base_name.str() + ".png");
  png.write(scanner.n_z_pixels, scanner.n_y_pixels, output);
}

void progress(int iteration, void* ptr, bool finished) {
  Context* context = static_cast<Context*>(ptr);
  context->progress(iteration, finished);
}

// wraps progress and output into abstract context ptr and run CUDA code
void run_reconstruction(Scanner<float, short>& scanner,
                        std::vector<Response<float>>& responses,
                        int n_iteration_blocks,
                        int n_iterations_per_block,
                        int device,
                        int n_blocks,
                        int n_threads_per_block,
                        bool verbose,
                        util::progress& progress,
                        std::string output_base_name,
                        bool text_output,
                        int n_file_digits) {

  Context context(progress, output_base_name, text_output, n_file_digits);
  run_reconstruction(scanner,
                     responses.data(),
                     responses.size(),
                     n_iteration_blocks,
                     n_iterations_per_block,
                     GPU::output,
                     GPU::progress,
                     &context,
                     device,
                     n_blocks,
                     n_threads_per_block,
                     verbose);
}

// convert double (DP) responses into float (SP) responses and run SP kernel
void run_reconstruction(Scanner<float, short>& scanner,
                        std::vector<Response<double>>& responses,
                        int n_iteration_blocks,
                        int n_iterations_per_block,
                        int device,
                        int n_blocks,
                        int n_threads_per_block,
                        bool verbose,
                        util::progress& progress,
                        std::string output_file_name,
                        bool text_output,
                        int n_file_digits) {
  std::vector<Response<float>> sp_responses;
  for (auto& response : responses) {
    Response<float> sp_response(response.z_u, response.z_d, response.dl);
    sp_responses.push_back(sp_response);
  }
  run_reconstruction(scanner,
                     sp_responses,
                     n_iteration_blocks,
                     n_iterations_per_block,
                     device,
                     n_blocks,
                     n_threads_per_block,
                     verbose,
                     progress,
                     output_file_name,
                     text_output,
                     n_file_digits);
}

}  // GPU
}  // Strip
}  // PET2D
