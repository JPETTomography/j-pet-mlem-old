#pragma once

#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>

#include "event.h"
#include "kernel.h"
#include "scanner.h"

#if _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#endif

#include "reconstruction_stats.h"

#define BB_UPDATE 1

namespace PET2D {
namespace Strip {

/// 2D strip PET reconstruction
template <typename FType = double, typename KernelType = Kernel<FType>>
class Reconstruction {
 public:
  using F = FType;
  using Kernel = KernelType;
  using Scanner = Strip::Scanner<FType, short>;
  using Pixel = typename Scanner::Pixel;
  using Point = typename Scanner::Point;
  using Vector = typename Scanner::Vector;

  Scanner scanner;

 private:
  const int n_threads;
  std::vector<F> rho;
  std::vector<F> acc_log;
  std::vector<std::vector<F>> thread_rhos;
  std::vector<F> sensitivity;
  Kernel kernel;
  ReconstructionStats<size_t> stats_;

 public:
  const ReconstructionStats<size_t>& stats;
  std::vector<Event<F>> events;

  Reconstruction(const Scanner& scanner)
      : scanner(scanner),
        n_threads(omp_get_max_threads()),
        rho(scanner.total_n_pixels, 100),
        thread_rhos(n_threads),
        sensitivity(scanner.total_n_pixels),
        kernel(scanner.sigma_z, scanner.sigma_dl),
        stats_(n_threads),
        stats(stats_) {

    for (int y = 0; y < scanner.n_y_pixels; ++y) {
      for (int z = 0; z < scanner.n_z_pixels; ++z) {

        sensitivity[y * scanner.n_z_pixels + z] =
            scanner.pixel_sensitivity(Pixel(z, y));
      }
    }
  }

  Reconstruction(F R_distance,
                 F scintilator_length,
                 int n_y_pixels,
                 int n_z_pixels,
                 F pixel_height,
                 F pixel_width,
                 F sigma_z,
                 F sigma_dl)
      : Reconstruction(Scanner(R_distance,
                               scintilator_length,
                               n_y_pixels,
                               n_z_pixels,
                               pixel_height,
                               pixel_width,
                               sigma_z,
                               sigma_dl)) {}

  Reconstruction(F R_distance,
                 F scintilator_length,
                 int n_pixels,
                 F pixel_size,
                 F sigma_z,
                 F sigma_dl)
      : Reconstruction(R_distance,
                       scintilator_length,
                       n_pixels,
                       n_pixels,
                       pixel_size,
                       pixel_size,
                       sigma_z,
                       sigma_dl) {}

  /// Performs n_iterations of the list mode MLEM algorithm
  template <typename ProgressCallback>
  void operator()(ProgressCallback progress,   ///< progress callback
                  int n_iterations,            ///< iterations to perform
                  int n_iterations_so_far = 0  ///< iterations so far
                  ) {

    stats_.fill();

    for (int iteration = 0; iteration < n_iterations; ++iteration) {
      progress(iteration + n_iterations_so_far);
      int n_events = events.size();

      for (auto& thread_rho : thread_rhos) {
        thread_rho.assign(scanner.total_n_pixels, 0);
      }

#if _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
      for (int e = 0; e < n_events; ++e) {
        int thread = omp_get_thread_num();
        stats_.n_events_processed_by(thread, 1);
#if BB_UPDATE
        auto event = events[e];
        F tan, y, z;
        event.transform(scanner.radius, tan, y, z);

        bb_update(Point(z, y), y, tan, thread_rhos[thread]);
#else
        F y = events[e].z_u;
        F z = events[e].z_d;
        simple_update(Point(z, y), y, z, thread_rhos[thread]);
#endif
      }

      rho.assign(scanner.total_n_pixels, 0);

      for (int thread = 0; thread < n_threads; ++thread) {
        for (int i = 0; i < scanner.total_n_pixels; ++i) {
          rho[i] += thread_rhos[thread][i];
        }
      }

      progress(iteration + n_iterations_so_far, true);
    }
    stats_.collect();
  }

  template <typename StreamType> Reconstruction& operator<<(StreamType& in) {
    int i = 0;
    while (!in.eof()) {

      F z_u, z_d, dl;
      in >> z_u >> z_d >> dl;
      // std::cout << i << " "<< z_u << " " << z_d << " " << dl << "\n";
      Event<F> temp_event(z_u, z_d, dl);
      events.push_back(temp_event);
      i++;
    }
    return *this;
  }

  template <typename StreamType> StreamType& operator>>(StreamType& out) {
    return out << rho;
  }

  template <class FileWriter>
  void output_bitmap(FileWriter& fw, bool output_sensitivity = false) {
    fw.template write_header<>(scanner.n_z_pixels, scanner.n_y_pixels);

    auto& output = output_sensitivity ? sensitivity : rho;
    F output_max = 0;
    for (auto& v : output) {
      output_max = std::max(output_max, v);
    }

    auto output_gain =
        static_cast<double>(std::numeric_limits<uint8_t>::max()) / output_max;

    uint8_t* row = (uint8_t*)alloca(scanner.n_z_pixels);
    for (int y = 0; y < scanner.n_y_pixels; ++y) {
      for (auto x = 0; x < scanner.n_z_pixels; ++x) {
        row[x] = std::numeric_limits<uint8_t>::max() -
                 output_gain * output[y * scanner.n_z_pixels + x];
      }
      fw.write_row(row);
    }
  }

 private:
  int n_pixels_in_line(F length, F pixel_size) const {
    return static_cast<int>(length / pixel_size + F(0.5));
  }

  void bb_update(Point ellipse_center, F y, F tan, std::vector<F>& output_rho) {
    bool use_sensitivity = false;
    F sec, A, B, C, bb_y, bb_z;
    kernel.ellipse_bb(tan, sec, A, B, C, bb_y, bb_z);

    Pixel center_pixel = scanner.pixel_at(ellipse_center);

    const int bb_half_width = n_pixels_in_line(bb_z, scanner.pixel_width);
    const int bb_half_height = n_pixels_in_line(bb_y, scanner.pixel_height);

    int thread = omp_get_thread_num();
    stats_.bb_width_sum_by(thread, 2 * bb_half_width);
    stats_.bb_height_sum_by(thread, 2 * bb_half_height);
    stats_.bb_width2_sum_by(thread, 4 * bb_half_width * bb_half_width);
    stats_.bb_height2_sum_by(thread, 4 * bb_half_height * bb_half_height);
    stats_.bb_width_height_sum_by(thread, 4 * bb_half_width * bb_half_height);

    Pixel top_left(center_pixel.x - bb_half_width,
                   center_pixel.y - bb_half_height);
    Pixel bottom_right(center_pixel.x + bb_half_width,
                       center_pixel.y + bb_half_height);
    Pixel scanner_top_left(0, 0);
    Pixel scanner_bottom_right(scanner.n_z_pixels - 1, scanner.n_y_pixels - 1);

    // check boundary conditions
    top_left.clamp(scanner_top_left, scanner_bottom_right);
    bottom_right.clamp(scanner_top_left, scanner_bottom_right);

    const int bb_size = 4 * bb_half_width * bb_half_height;
    F* ellipse_kernel_mul_rho = (F*)alloca(bb_size * sizeof(F));
    Pixel* ellipse_pixels = (Pixel*)alloca(bb_size * sizeof(Pixel));
    int n_ellipse_pixels = 0;

    F denominator = 0;

    for (int iy = top_left.y; iy < bottom_right.y; ++iy) {
      for (int iz = top_left.x; iz < bottom_right.x; ++iz) {
        stats_.n_pixels_processed_by();
        Pixel pixel(iz, iy);
        Point point = scanner.pixel_center(pixel);

        if (kernel.in_ellipse(A, B, C, ellipse_center, point)) {
          Vector r = point - ellipse_center;

          int i = pixel.y * scanner.n_z_pixels + pixel.x;

          F pixel_sensitivity = use_sensitivity ? sensitivity[i] : 1;
          stats_.n_kernel_calls_by();
          F event_kernel = kernel(y, tan, sec, scanner.radius, r);
          F event_kernel_mul_rho = event_kernel * rho[i];
          denominator += event_kernel_mul_rho * pixel_sensitivity;
          ellipse_pixels[n_ellipse_pixels] = pixel;
          ellipse_kernel_mul_rho[n_ellipse_pixels] = event_kernel_mul_rho;
          ++n_ellipse_pixels;
        }
      }
    }

    F inv_denominator = 1 / denominator;

    for (int p = 0; p < n_ellipse_pixels; ++p) {
      auto pixel = ellipse_pixels[p];
      auto pixel_kernel = ellipse_kernel_mul_rho[p];
      int i = pixel.y * scanner.n_z_pixels + pixel.x;
      output_rho[i] += pixel_kernel * inv_denominator;
    }
  }

  void simple_update(Point ellipse_center,
                     F y,
                     F z,
                     std::vector<F>& output_rho) {

    Pixel center_pixel = scanner.pixel_at(ellipse_center);

    int y_line = 3 * scanner.sigma_z / scanner.pixel_width;
    int z_line = 3 * scanner.sigma_dl / scanner.pixel_height;

    const Pixel tl(center_pixel.x - z_line, center_pixel.y - y_line);
    const Pixel br(center_pixel.x + z_line, center_pixel.y + y_line);

    const int bb_size = 4 * z_line * y_line;
    F* ellipse_kernel_mul_rho = (F*)alloca(bb_size * sizeof(F));
    int n_ellipse_pixels = 0;

    F denominator = 0;

    for (int iy = tl.y; iy < br.y; ++iy) {
      for (int iz = tl.x; iz < br.x; ++iz) {
        stats_.n_pixels_processed_by();
        Pixel pixel(iz, iy);
        Point point = scanner.pixel_center(pixel);

        int i = pixel.y * scanner.n_z_pixels + pixel.x;
        stats_.n_kernel_calls_by();
        F event_kernel =
            kernel.test(y, z, point, scanner.sigma_z, scanner.sigma_dl);
        F event_kernel_mul_rho = event_kernel * rho[i];
        ellipse_kernel_mul_rho[n_ellipse_pixels++] = event_kernel_mul_rho;
      }
    }

    F inv_denominator = 1 / denominator;

    for (int iy = tl.y; iy < br.y; ++iy) {
      for (int iz = tl.x; iz < br.x; ++iz) {
        Pixel pixel(iz, iy);
        int i = pixel.y * scanner.n_z_pixels + pixel.x;
        int ik = (pixel.y - tl.y) * z_line * 2 + pixel.x - tl.x;
        output_rho[i] += ellipse_kernel_mul_rho[ik] * rho[i] * inv_denominator;
      }
    }
  }
};
}  // Strip
}  // PET2D
