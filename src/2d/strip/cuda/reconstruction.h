#pragma once

#include <vector>
#include <sstream>
#include <iomanip>

#include "util/png_writer.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"
#include "util/progress.h"

#include "../event.h"
#include "../strip_detector.h"

// this is forward declaration of CUDA side entry-point function
template <typename F>
void run_gpu_reconstruction(StripDetector<F>& detector,
                            Event<F>* events,
                            int n_events,
                            int n_iteration_blocks,
                            int n_iterations_in_block,
                            void (*output_callback)(StripDetector<F>& detector,
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

namespace GPU {

struct Context {
  Context(Progress& progress, std::string& output_file_name)
      : progress(progress), output_file_name(output_file_name) {}
  Progress& progress;
  std::string& output_file_name;
};

void output(StripDetector<float>& detector,
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
    base_name << std::setw(3) << std::setfill('0')  //
              << iteration + 1 << std::setw(0);     // 001
  } else {
    base_name << "sensitivity";
  }

  obstream bin(base_name.str() + ".bin");
  bin.write(output, detector.total_n_pixels);

  png_writer png(base_name.str() + ".png");
  png.write_header<>(detector.n_z_pixels, detector.n_y_pixels);

  float output_max = 0;
  for (int i = 0; i < detector.total_n_pixels; ++i) {
    output_max = std::max(output_max, output[i]);
  }

  auto output_gain =
      static_cast<double>(std::numeric_limits<uint8_t>::max()) / output_max;

  uint8_t* row = (uint8_t*)alloca(detector.n_z_pixels);
  for (int y = 0; y < detector.n_y_pixels; ++y) {
    for (auto x = 0; x < detector.n_z_pixels; ++x) {
      row[x] = std::numeric_limits<uint8_t>::max() -
               output_gain * output[y * detector.n_z_pixels + x];
    }
    png.write_row(row);
  }
}

void progress(int iteration, void* ptr, bool finished) {
  Context* context = static_cast<Context*>(ptr);
  context->progress(iteration, finished);
}
}

// wraps progress and output into abstract context ptr and run CUDA code
void run_gpu_reconstruction(StripDetector<float>& detector,
                            std::vector<Event<float>>& events,
                            int n_iteration_blocks,
                            int n_iterations_per_block,
                            int device,
                            int n_blocks,
                            int n_threads_per_block,
                            bool verbose,
                            Progress& progress,
                            std::string output) {

  GPU::Context context(progress, output);
  run_gpu_reconstruction(detector,
                         events.data(),
                         events.size(),
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

// convert double (DP) events into float (SP) events and run SP kernel
void run_gpu_reconstruction(StripDetector<float>& detector,
                            std::vector<Event<double>>& events,
                            int n_iteration_blocks,
                            int n_iterations_per_block,
                            int device,
                            int n_blocks,
                            int n_threads_per_block,
                            bool verbose,
                            Progress& progress,
                            std::string output_file_name) {
  std::vector<Event<float>> sp_event_list;
  for (auto& event : events) {
    Event<float> sp_event(event.z_u, event.z_d, event.dl);
    sp_event_list.push_back(sp_event);
  }
  run_gpu_reconstruction(detector,
                         sp_event_list,
                         n_iteration_blocks,
                         n_iterations_per_block,
                         device,
                         n_blocks,
                         n_threads_per_block,
                         verbose,
                         progress,
                         output_file_name);
}
