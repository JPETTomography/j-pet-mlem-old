#pragma once

#include <vector>

#include "util/png_writer.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"

#include "../event.h"
#include "../strip_detector.h"

#include "config.h"

template <typename F>
void run_reconstruction_kernel(StripDetector<F>& detector,
                               Event<F>* events,
                               int events_size,
                               int iteration_chunk,
                               F* image_output,
                               int device,
                               int n_blocks,
                               int n_threads_per_block);

void run_gpu_reconstruction(StripDetector<float>& detector,
                            std::vector<Event<float>>& events,
                            int device,
                            int n_blocks,
                            int n_threads_per_block) {

  std::vector<std::vector<float>> image_output;
  image_output.assign(detector.n_y_pixels,
                      std::vector<float>(detector.n_y_pixels, float(0)));

  std::vector<float> gpu_output_image;
  gpu_output_image.resize(n_blocks * detector.total_n_pixels);
  image_output.assign(n_blocks * detector.n_y_pixels,
                      std::vector<float>(detector.n_y_pixels, float(0)));

  run_reconstruction_kernel(detector,
                            events.data(),
                            events.size(),
                            n_blocks,
                            gpu_output_image.data(),
                            device,
                            n_blocks,
                            n_threads_per_block);

  for (int iteration = 0; iteration < n_blocks; ++iteration) {
    for (int i = 0; i < detector.n_y_pixels; ++i) {
      for (int j = 0; j < detector.n_z_pixels; ++j) {
        image_output[i][j] =
            gpu_output_image[iteration * detector.total_n_pixels +
                             (i * detector.n_z_pixels + j)];
      }
    }

    png_writer png("gpu_rec_i_" + std::to_string(iteration + 1) + ".png");
    png.write_header<>(detector.n_y_pixels, detector.n_z_pixels);

    float output_max = 0.0;
    for (auto& col : image_output) {
      for (auto& row : col) {
        output_max = std::max(output_max, row);
      }
    }

    std::ofstream data_output("pixels_output_i_" +
                              std::to_string(iteration + 1) + ".txt");

    for (int x = 0; x < detector.n_y_pixels; ++x) {
      for (int y = 0; y < detector.n_z_pixels; ++y) {
        data_output << x << " " << y << " " << image_output[x][y] << std::endl;
      }
    }

    auto output_gain =
        static_cast<double>(std::numeric_limits<uint8_t>::max()) / output_max;

    for (int y = 0; y < detector.n_y_pixels; ++y) {
      uint8_t row[detector.n_z_pixels];
      for (auto x = 0; x < detector.n_z_pixels; ++x) {
        row[x] = std::numeric_limits<uint8_t>::max() -
                 output_gain * image_output[y][x];
      }
      png.write_row(row);
    }
  }
}

void run_gpu_reconstruction(StripDetector<float>& detector,
                            std::vector<Event<double>>& events,
                            int device,
                            int n_blocks,
                            int n_threads_per_block) {
  std::vector<Event<float>> sp_event_list;
  for (auto& event : events) {
    Event<float> sp_event(event.z_u, event.z_d, event.dl);
    sp_event_list.push_back(sp_event);
  }
  run_gpu_reconstruction(
      detector, sp_event_list, device, n_blocks, n_threads_per_block);
}
