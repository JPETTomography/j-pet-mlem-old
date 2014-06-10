#include <cmath>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <unordered_map>
#include "event.h"

#include "util/png_writer.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"

void gpu_reconstruction_strip_2d(gpu_config::GPU_parameters cfg,
                                 event<float>* event_list,
                                 int event_size,
                                 int iteration_chunk,
                                 float* image_output,
                                 int off);

void execute_kernel_reconstruction(gpu_config::GPU_parameters cfg,
                                   event<float>* event_list,
                                   int event_size,
                                   int warp_offset,
                                   int n_blocks) {

  std::vector<std::vector<float>> image_output;
  image_output.assign(cfg.n_pixels, std::vector<float>(cfg.n_pixels, float(0)));

  std::vector<float> gpu_output_image;
  gpu_output_image.resize(n_blocks * cfg.n_pixels * cfg.n_pixels);
  image_output.assign(n_blocks * cfg.n_pixels,
                      std::vector<float>(cfg.n_pixels, float(0)));

  gpu_reconstruction_strip_2d(cfg,
                              event_list,
                              event_size,
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

    std::string file = std::string("gpu_rec_i_");

    file.append(std::to_string(iteration + 1));
    file.append(".png");

    png_writer png(file);
    png.write_header<>(cfg.n_pixels, cfg.n_pixels);

    float output_max = 0.0;
    for (auto& col : image_output) {
      for (auto& row : col) {
        output_max = std::max(output_max, row);
      }
    }

    std::ofstream data_output;
    data_output.open("pixels_output.txt");
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
