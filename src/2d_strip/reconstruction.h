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

#define BB_UPDATE 1

template <typename FType = double> class Reconstruction {
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
  std::vector<F> thread_rho;
  std::vector<F> sensitivity;

  Kernel<F> kernel;

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
        thread_rho(detector.total_n_pixels * n_threads),
        sensitivity(detector.total_n_pixels) {

    for (int y = 0; y < detector.n_y_pixels; ++y) {
      for (int z = 0; z < detector.n_z_pixels; ++z) {
        Point point = detector.pixel_center(Pixel(z, y));
        sensitivity[y * detector.n_z_pixels + z] = detector.sensitivity(point);
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

    for (int iteration = 0; iteration < n_iterations; ++iteration) {
      thread_rho.assign(detector.total_n_pixels * n_threads, 0);
      progress(iteration + n_iterations_so_far);
      int n_events = events.size();

#if _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
      for (int e = 0; e < n_events; ++e) {
        int thread = omp_get_thread_num();

#if BB_UPDATE
        auto event = events[e];
        F tan, y, z;
        event.transform(detector.radius, tan, y, z);
        F angle = std::atan(tan);

        bb_update(Point(z, y), angle, y, tan, thread);
#else
        F y = events[id].z_u;
        F z = events[id].z_d;
        simple_update(Point(z, y), y, z, tid, id);
#endif
      }

      rho.assign(detector.total_n_pixels, 0);

      for (int thread = 0; thread < n_threads; ++thread) {
        for (int i = 0; i < detector.total_n_pixels; ++i) {
          rho[i] += thread_rho[thread * detector.total_n_pixels + i];
        }
      }
    }
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

  void bb_update(Point ellipse_center, F angle, F y, F tan, int thread) {

    F sec, sec_sq, A, B, C, bb_y, bb_z;
    detector.ellipse_bb(angle, tan, sec, sec_sq, A, B, C, bb_y, bb_z);

    Pixel center_pixel = detector.pixel_location(ellipse_center);

    const int bb_half_width = n_pixels_in_line(bb_z, detector.pixel_width);
    const int bb_half_height = n_pixels_in_line(bb_y, detector.pixel_height);
    const Pixel tl(center_pixel.x - bb_half_width,
                   center_pixel.y - bb_half_height);
    const Pixel br(center_pixel.x + bb_half_width,
                   center_pixel.y + bb_half_height);

    std::vector<std::pair<Pixel, F>> ellipse_kernel_mul_rho;
    ellipse_kernel_mul_rho.reserve(2000);

    F acc = 0;

    for (int iy = tl.y; iy < br.y; ++iy) {
      for (int iz = tl.x; iz < br.x; ++iz) {
        Pixel pixel(iz, iy);
        Point point = detector.pixel_center(pixel);

        if (detector.in_ellipse(A, B, C, ellipse_center, point)) {
          point -= ellipse_center;

          int i = pixel.y * detector.n_y_pixels + pixel.x;

          F pixel_sensitivity = sensitivity[i];
          F event_kernel = kernel(y,
                                  tan,
                                  sec,
                                  sec_sq,
                                  detector.radius,
                                  point,
                                  detector.inv_cor_mat_diag,
                                  sqrt_det_cor_mat) /
                           pixel_sensitivity;
          F event_kernel_mul_rho = event_kernel * rho[i];
          ellipse_kernel_mul_rho.push_back(
              std::make_pair(pixel, event_kernel_mul_rho));
          acc += event_kernel_mul_rho * sensitivity[i];
        }
      }
    }

    F inv_acc = 1 / acc;

    for (auto& e : ellipse_kernel_mul_rho) {
      auto pixel = e.first;
      int i = pixel.y * detector.n_y_pixels + pixel.x;
      thread_rho[thread * detector.total_n_pixels + i] += e.second * inv_acc;
    }
  }

  void simple_update(Point ellipse_center, F y, F z, int thread) {

    Pixel center_pixel = detector.pixel_location(ellipse_center);

    int y_line = 3 * (detector.sigma_z / detector.pixel_width);
    int z_line = 3 * (detector.sigma_dl / detector.pixel_height);

    std::vector<std::pair<Pixel, F>> ellipse_kernels;
    ellipse_kernels.reserve(2000);

    F acc = 0;

    for (int iy = center_pixel.y - y_line; iy < center_pixel.y + y_line; ++iy) {
      for (int iz = center_pixel.x - z_line; iz < center_pixel.x + z_line;
           ++iz) {
        int i = iy * detector.n_z_pixels + iz;
        Pixel pixel(iz, iy);
        Point point = detector.pixel_center(pixel);

        F event_kernel = kernel.test_kernel(
            y, z, point, detector.sigma_z, detector.sigma_dl);

        ellipse_kernels.push_back(std::make_pair(pixel, event_kernel));
        acc += event_kernel * rho[i];
      }
    }

    for (auto& e : ellipse_kernels) {
      auto pixel = e.first;
      int i = pixel.y * detector.n_z_pixels + pixel.x;
      thread_rho[thread * detector.total_n_pixels + i] +=
          e.second * rho[i] / acc;
    }
  }
};
