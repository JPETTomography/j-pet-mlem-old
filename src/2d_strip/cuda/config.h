#pragma once

#define WARP_SIZE 32
//#define WARP_GRANULARITY
#define EVENT_GRANULARITY
#define NORMAL_PHANTOM 0

namespace CUDA {
struct Config {
  float R_distance, Scentilator_length, pixel_size;
  int n_pixels;
  float sigma;
  float dl;
  int number_of_blocks;
  int number_of_threads_per_block;
  int number_of_events;
  float inv_pow_sigma_dl;
  float inv_pow_sigma_z;
  float grid_size_y_;
  float grid_size_z_;
};
}
