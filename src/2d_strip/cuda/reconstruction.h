#pragma once

#include <vector>
#include <sstream>
#include <iomanip>

#include "util/png_writer.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"
#include "util/util.h"

#include "../event.h"
#include "../strip_detector.h"

#include "config.h"

template <typename F>
void run_reconstruction_kernel(
    StripDetector<F>& detector,
    Event<F>* events,
    int n_events,
    int n_iteration_blocks,
    int n_iterations_in_block,
    void (*output_callback)(StripDetector<F>& detector,
                            int iteration,
                            F* image,
                            void* context),
    void (*progress_callback)(int iteration, void* context),
    void* context,
    int device,
    int n_blocks,
    int n_threads_per_block);

namespace GPU {
struct Context {
  Context(Progress& progress, std::string& output)
      : progress(progress), output(output) {}
  Progress& progress;
  std::string& output;
};

void output(StripDetector<float>& detector,
            int iteration,
            float* output,
            void* ptr) {
  Context* context = static_cast<Context*>(ptr);
  std::stringstream fn;
  fn << context->output << "_";
  if (iteration >= 0) {
    fn << std::setw(3) << std::setfill('0') << iteration << std::setw(0)
       << ".png";
  } else {
    fn << "sensitivity.png";
  }
  png_writer png(fn.str());
  png.write_header<>(detector.n_y_pixels, detector.n_z_pixels);

  float output_max = 0;
  for (int i = 0; i < detector.total_n_pixels; ++i) {
    output_max = std::max(output_max, output[i]);
  }

  auto output_gain =
      static_cast<double>(std::numeric_limits<uint8_t>::max()) / output_max;

  for (int y = 0; y < detector.n_y_pixels; ++y) {
    uint8_t row[detector.n_z_pixels];
    for (auto x = 0; x < detector.n_z_pixels; ++x) {
      row[x] = std::numeric_limits<uint8_t>::max() -
               output_gain * output[y * detector.n_z_pixels + x];
    }
    png.write_row(row);
  }
}

void progress(int iteration, void* ptr) {
  Context* context = static_cast<Context*>(ptr);
  context->progress(iteration);
}
}

void run_gpu_reconstruction(StripDetector<float>& detector,
                            std::vector<Event<float>>& events,
                            int n_iteration_blocks,
                            int n_iterations_per_block,
                            int device,
                            int n_blocks,
                            int n_threads_per_block,
                            Progress& progress,
                            std::string output) {

  GPU::Context context(progress, output);
  run_reconstruction_kernel(detector,
                            events.data(),
                            events.size(),
                            n_iteration_blocks,
                            n_iterations_per_block,
                            GPU::output,
                            GPU::progress,
                            &context,
                            device,
                            n_blocks,
                            n_threads_per_block);
}

void run_gpu_reconstruction(StripDetector<float>& detector,
                            std::vector<Event<double>>& events,
                            int n_iteration_blocks,
                            int n_iterations_per_block,
                            int device,
                            int n_blocks,
                            int n_threads_per_block,
                            Progress& progress,
                            std::string output) {
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
                         progress,
                         output);
}
