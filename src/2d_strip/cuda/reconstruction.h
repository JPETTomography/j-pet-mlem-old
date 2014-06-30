#pragma once

#include <cmath>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <unordered_map>

#include "util/png_writer.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"

#include "../event.h"

#include "config.h"

void gpu_reconstruction_strip_2d(CUDA::Config cfg,
                                 Event<float>* event_list,
                                 int event_size,
                                 int iteration_chunk,
                                 float* image_output,
                                 int off);

void execute_kernel_reconstruction(CUDA::Config cfg,
                                   std::vector<Event<float>>& event_list,
                                   int warp_offset,
                                   int n_blocks) {

  std::vector<std::vector<float>> image_output;
  image_output.assign(cfg.n_pixels, std::vector<float>(cfg.n_pixels, float(0)));

  std::vector<float> gpu_output_image;
  gpu_output_image.resize(n_blocks * cfg.n_pixels * cfg.n_pixels);
  image_output.assign(n_blocks * cfg.n_pixels,
                      std::vector<float>(cfg.n_pixels, float(0)));

  gpu_reconstruction_strip_2d(cfg,
                              event_list.data(),
                              event_list.size(),
                              n_blocks,
                              gpu_output_image.data(),
                              warp_offset);

  for (int iteration = 0; iteration < n_blocks; ++iteration) {

    int mem_offset = cfg.n_pixels * cfg.n_pixels;

    for (int i = 0; i < cfg.n_pixels; ++i) {
      for (int j = 0; j < cfg.n_pixels; ++j) {

        image_output[i][j] =
            gpu_output_image[iteration * mem_offset + (i * cfg.n_pixels + j)];
      }
    }

    png_writer png("gpu_rec_i_" + std::to_string(iteration + 1) + ".png");
    png.write_header<>(cfg.n_pixels, cfg.n_pixels);

    float output_max = 0.0;
    for (auto& col : image_output) {
      for (auto& row : col) {
        output_max = std::max(output_max, row);
      }
    }

    std::ofstream data_output("pixels_output_i_" +
                              std::to_string(iteration + 1) + ".txt");

    for (int x = 0; x < cfg.n_pixels; ++x) {
      for (int y = 0; y < cfg.n_pixels; ++y) {

        data_output << x << " " << y << " " << image_output[x][y] << std::endl;
      }
    }

    auto output_gain =
        static_cast<double>(std::numeric_limits<uint8_t>::max()) / output_max;

    for (int y = 0; y < cfg.n_pixels; ++y) {
      uint8_t row[cfg.n_pixels];
      for (auto x = 0; x < cfg.n_pixels; ++x) {
        row[x] = std::numeric_limits<uint8_t>::max() -
                 output_gain * image_output[y][x];
      }
      png.write_row(row);
    }
  }
}

void execute_kernel_reconstruction(CUDA::Config cfg,
                                   std::vector<Event<double>>& event_list,
                                   int warp_offset,
                                   int n_blocks) {
  std::vector<Event<float>> sp_event_list;
  for (auto& event : event_list) {
    Event<float> sp_event(event.z_u, event.z_d, event.dl);
    sp_event_list.push_back(sp_event);
  }
  execute_kernel_reconstruction(cfg, sp_event_list, warp_offset, n_blocks);
}
