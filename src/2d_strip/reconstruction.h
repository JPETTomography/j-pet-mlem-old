#pragma once

#include <cmath>
#include <vector>
#include <algorithm>

#include "event.h"
#include "kernel.h"
#include "strip_detector.h"

#if _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#endif

#include "reconstruction_stats.h"

#define BB_UPDATE 1


#define READER(T, name) \
  T name() { return stats_.total_##name##_; }

template <typename FType = double, template <typename Ft> class K = Kernel>
class Reconstruction {
 public:
  typedef FType F;
  typedef StripDetector<FType> Detector;
  typedef typename Detector::Pixel Pixel;
  typedef typename Detector::Point Point;

  Detector detector;
  F sqrt_det_cor_mat;

 private:
  const int n_threads;
  std::vector<Event<F>> events;
  std::vector<F> rho;
  std::vector<F> acc_log;
  std::vector<std::vector<F>> thread_rhos;
  std::vector<F> sensitivity;
  std::vector<F> inv_sensitivity;

  K<F> kernel;

  Stats<size_t> stats_;

 public:
  Reconstruction(F R_distance,
                 F scintilator_length,
                 int n_y_pixels,
                 int n_z_pixels,
                 F pixel_height,
                 F pixel_width,
                 F sigma_z,
                 F sigma_dl)
      : detector(R_distance,
                 scintilator_length,
                 n_y_pixels,
                 n_z_pixels,
                 pixel_height,
                 pixel_width,
                 sigma_z,
                 sigma_dl),
        sqrt_det_cor_mat(detector.sqrt_det_cor_mat()),
        n_threads(omp_get_max_threads()),
        rho(detector.total_n_pixels, 100),
        thread_rhos(n_threads),
        sensitivity(detector.total_n_pixels),
        inv_sensitivity(detector.total_n_pixels),
        stats_(n_threads) {

    for (int y = 0; y < detector.n_y_pixels; ++y) {
      for (int z = 0; z < detector.n_z_pixels; ++z) {
#ifdef USE_SENSITIVITY
        Point point = detector.pixel_center(Pixel(z, y));
        F pixel_sensitivity = detector.sensitivity(point)/3;
        Point ur(pixel_width / 2, pixel_height / 2);
        pixel_sensitivity += detector.sensitivity(point + ur)/6;
        pixel_sensitivity += detector.sensitivity(point - ur)/6;
        Point ul(-pixel_width / 2, pixel_height / 2);
        pixel_sensitivity += detector.sensitivity(point + ul)/6;
        pixel_sensitivity += detector.sensitivity(point - ul)/6;
        sensitivity[y * detector.n_z_pixels + z] = pixel_sensitivity;
        inv_sensitivity[y * detector.n_z_pixels + z] = pixel_sensitivity;

#else
        sensitivity[y * detector.n_z_pixels + z] = 1;
        inv_sensitivity[y * detector.n_z_pixels + z] = 1;
#endif
      }
    }
  }

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
  void operator()(ProgressCallback progress,
                  int n_iterations,
                  int n_iterations_so_far = 0) {

    stats_.fill(0);

    for (int iteration = 0; iteration < n_iterations; ++iteration) {
      for (auto& rho : thread_rhos) {
        rho.assign(detector.total_n_pixels, 0);
      }
      progress(iteration + n_iterations_so_far);
      int n_events = events.size();

#if _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
      for (int e = 0; e < n_events; ++e) {
        int thread = omp_get_thread_num();
        stats_.n_events_processed_[thread]++;
#if BB_UPDATE
        auto event = events[e];
        F tan, y, z;
        event.transform(detector.radius, tan, y, z);
        F angle = std::atan(tan);

        bb_update(Point(z, y), angle, y, tan, thread_rhos[thread]);
#else
        F y = events[e].z_u;
        F z = events[e].z_d;
        simple_update(Point(z, y), y, z, thread_rhos[thread]);
#endif
      }

      rho.assign(detector.total_n_pixels, 0);

      for (int thread = 0; thread < n_threads; ++thread) {
        for (int i = 0; i < detector.total_n_pixels; ++i) {
          rho[i] += thread_rhos[thread][i];
        }
      }
    }
    stats_.collect();
  }

  // accessor for CUDA compatibility
  std::vector<Event<F>>& get_event_list() { return events; }

  template <typename StreamType> Reconstruction& operator<<(StreamType& in) {
    int size;
    in >> size;

    for (int it = 0; it < size; ++it) {
      F z_u, z_d, dl;
      in >> z_u >> z_d >> dl;
      Event<F> temp_event(z_u, z_d, dl);
      events.push_back(temp_event);
    }
    return *this;
  }

  template <class FileWriter>
  void output_bitmap(FileWriter& fw, bool output_sensitivity = false) {
    fw.template write_header<>(detector.n_z_pixels, detector.n_y_pixels);

    auto& output = output_sensitivity ? sensitivity : rho;
    F output_max = 0;
    for (auto& v : output) {
      output_max = std::max(output_max, v);
    }

    auto output_gain =
        static_cast<double>(std::numeric_limits<uint8_t>::max()) / output_max;

    uint8_t* row = (uint8_t*)alloca(detector.n_z_pixels);
    for (int y = 0; y < detector.n_y_pixels; ++y) {
      for (auto x = 0; x < detector.n_z_pixels; ++x) {
        row[x] = std::numeric_limits<uint8_t>::max() -
                 output_gain * output[y * detector.n_z_pixels + x];
      }
      fw.write_row(row);
    }
  }

 private:
  int n_pixels_in_line(F length, F pixel_size) const {
    return static_cast<int>((length + F(0.5)) / pixel_size);
  }

  void bb_update(Point ellipse_center,
                 F angle,
                 F y,
                 F tan,
                 std::vector<F>& output_rho) {

    F sec, sec_sq, A, B, C, bb_y, bb_z;
    detector.ellipse_bb(angle, tan, sec, sec_sq, A, B, C, bb_y, bb_z);

    Pixel center_pixel = detector.pixel_location(ellipse_center);

    const int bb_half_width = n_pixels_in_line(bb_z, detector.pixel_width);
    const int bb_half_height = n_pixels_in_line(bb_y, detector.pixel_height);

    int thread_number = omp_get_thread_num();
    stats_.bb_width_sum_[thread_number] += 2 * bb_half_width;
    stats_.bb_height_sum_[thread_number] += 2 * bb_half_height;
    stats_.bb_width2_sum_[thread_number] += 4 * bb_half_width * bb_half_width;
    stats_.bb_height2_sum_[thread_number] +=
        4 * bb_half_height * bb_half_height;
    stats_.bb_width_height_sum_[thread_number] +=
        4 * bb_half_width * bb_half_height;

    const Pixel tl(center_pixel.x - bb_half_width,
                   center_pixel.y - bb_half_height);
    const Pixel br(center_pixel.x + bb_half_width,
                   center_pixel.y + bb_half_height);

    const int bb_size = 4 * bb_half_width * bb_half_height;
    F* ellipse_kernel_mul_rho = (F*)alloca(bb_size * sizeof(F));
    Pixel* ellipse_pixels = (Pixel*)alloca(bb_size * sizeof(Pixel));
    int n_ellipse_pixels = 0;

    F acc = 0;

    for (int iy = tl.y; iy < br.y; ++iy) {
      for (int iz = tl.x; iz < br.x; ++iz) {
        stats_.n_pixels_processed_[omp_get_thread_num()]++;
        Pixel pixel(iz, iy);
        Point point = detector.pixel_center(pixel);

        if (detector.in_ellipse(A, B, C, ellipse_center, point)) {
          point -= ellipse_center;

          int i = pixel.y * detector.n_z_pixels + pixel.x;

          F pixel_sensitivity = sensitivity[i];
          stats_.n_kernel_calls_[omp_get_thread_num()]++;
          F event_kernel = kernel(y,
                                  tan,
                                  sec,
                                  sec_sq,
                                  detector.radius,
                                  point,
                                  detector.inv_cor_mat_diag,
                                  sqrt_det_cor_mat);
          F event_kernel_mul_rho = event_kernel * rho[i];
          acc += event_kernel_mul_rho*pixel_sensitivity;
          //event_kernel_mul_rho *= pixel_sensitivity;
          ellipse_pixels[n_ellipse_pixels] = pixel;
          ellipse_kernel_mul_rho[n_ellipse_pixels] = event_kernel_mul_rho;
          ++n_ellipse_pixels;
        }
      }
    }

    F inv_acc = 1 / acc;

    for (int p = 0; p < n_ellipse_pixels; ++p) {
      auto pixel = ellipse_pixels[p];
      auto pixel_kernel = ellipse_kernel_mul_rho[p];
      int i = pixel.y * detector.n_z_pixels + pixel.x;
      output_rho[i] += pixel_kernel * inv_acc;
    }
  }

  void simple_update(Point ellipse_center,
                     F y,
                     F z,
                     std::vector<F>& output_rho) {

    Pixel center_pixel = detector.pixel_location(ellipse_center);

    int y_line = 3 * detector.sigma_z / detector.pixel_width;
    int z_line = 3 * detector.sigma_dl / detector.pixel_height;

    const Pixel tl(center_pixel.x - z_line, center_pixel.y - y_line);
    const Pixel br(center_pixel.x + z_line, center_pixel.y + y_line);

    const int bb_size = 4 * z_line * y_line;
    F* ellipse_kernel_mul_rho = (F*)alloca(bb_size * sizeof(F));
    int n_ellipse_pixels = 0;

    F acc = 0;

    for (int iy = tl.y; iy < br.y; ++iy) {
      for (int iz = tl.x; iz < br.x; ++iz) {
        stats_.n_pixels_processed_[omp_get_thread_num()]++;
        Pixel pixel(iz, iy);
        Point point = detector.pixel_center(pixel);

        int i = pixel.y * detector.n_z_pixels + pixel.x;
        stats_.n_kernel_calls_[omp_get_thread_num()]++;
        F event_kernel =
            kernel.test(y, z, point, detector.sigma_z, detector.sigma_dl);
        F event_kernel_mul_rho = event_kernel * rho[i];
        ellipse_kernel_mul_rho[n_ellipse_pixels++] = event_kernel_mul_rho;
      }
    }

    F inv_acc = 1 / acc;

    for (int iy = tl.y; iy < br.y; ++iy) {
      for (int iz = tl.x; iz < br.x; ++iz) {
        Pixel pixel(iz, iy);
        int i = pixel.y * detector.n_z_pixels + pixel.x;
        int ik = (pixel.y - tl.y) * z_line * 2 + pixel.x - tl.x;
        output_rho[i] += ellipse_kernel_mul_rho[ik] * rho[i] * inv_acc;
      }
    }
  }

 public:
  READER(size_t, n_events_processed);
  READER(size_t, n_pixels_processed);
  READER(size_t, n_kernel_calls);
  READER(size_t, bb_width_sum);
  READER(size_t, bb_width2_sum);
  READER(size_t, bb_height_sum);
  READER(size_t, bb_height2_sum);
  READER(size_t, bb_width_height_sum);
};
