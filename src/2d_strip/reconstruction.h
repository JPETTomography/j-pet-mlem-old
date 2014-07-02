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

template <typename FType = double, typename DetectorType = StripDetector<FType>>
class Reconstruction {
 public:
  typedef FType F;
  typedef DetectorType Detector;
  typedef typename Detector::Pixel Pixel;
  typedef typename Detector::Point Point;

  F sqrt_det_cor_mat;
  Detector detector;

 private:
  std::vector<Event<F>> events;
  std::vector<std::vector<F>> rho;
  std::vector<std::vector<F>> rho_temp;
  std::vector<F> acc_log;
  std::vector<std::vector<F>> thread_rho;
  std::vector<std::vector<F>> sensitivity;

  Kernel<F> kernel;

  const int n_pixels;
  const int pixel_size;

 public:
  Reconstruction(F R_distance,
                 F scintilator_length,
                 int n_pixels,
                 F pixel_size,
                 F sigma_z,
                 F sigma_dl)
      : detector(R_distance,
                 scintilator_length,
                 n_pixels,
                 n_pixels,
                 pixel_size,
                 pixel_size,
                 sigma_z,
                 sigma_dl),
        n_pixels(n_pixels),
        pixel_size(pixel_size) {
    kernel = Kernel<F>();
    rho.assign(n_pixels, std::vector<F>(n_pixels, 100));
    rho_temp.assign(n_pixels, std::vector<F>(n_pixels, 10));

    sqrt_det_cor_mat = std::sqrt(detector.inv_cor_mat_diag[0] *  //
                                 detector.inv_cor_mat_diag[1] *  //
                                 detector.inv_cor_mat_diag[2]);

    sensitivity.assign(n_pixels, std::vector<F>(n_pixels));

    for (int y = 0; y < n_pixels; ++y) {
      for (int z = 0; z < n_pixels; ++z) {
        Point point = detector.pixel_center(y, z);
        sensitivity[y][z] = detector.sensitivity(point.x, point.y);
      }
    }
  }

  /// Performs n_iterations of the list mode MLEM algorithm
  template <typename ProgressCallback>
  void operator()(ProgressCallback progress,
                  int n_iterations,
                  int n_iterations_so_far = 0) {

    for (int i = 0; i < n_iterations; i++) {

      thread_rho.assign(omp_get_max_threads(),
                        std::vector<F>(n_pixels * n_pixels, 0));

      progress(i + n_iterations_so_far);

      int size = events.size();

#if _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
      for (int id = 0; id < size; ++id) {
        int tid = omp_get_thread_num();

#if BB_UPDATE
        auto event = events[id];
        F tan = event.tan(detector.radius);
        F y = event.y(tan);
        F z = event.z(y, tan);
        F angle = std::atan(tan);

        Point ellipse_center = Point(y, z);
        bb_update(ellipse_center, angle, y, tan, tid);
#else
        F y = events[id].z_u;
        F z = events[id].z_d;

        Point ellipse_center = Point(y, z);
        simple_update(ellipse_center, y, z, tid, id);
#endif
      }

      rho.assign(n_pixels, std::vector<F>(n_pixels, 0));

      for (int i = 0; i < n_pixels; ++i) {
        for (int j = 0; j < n_pixels; ++j) {
          for (int k = 0; k < omp_get_max_threads(); ++k) {
            rho[i][j] += thread_rho[k][i + j * n_pixels];
          }
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

  template <class FileWriter> void output_bitmap(FileWriter& fw) {
    fw.template write_header<>(n_pixels, n_pixels);

    F output_max = 0;
    for (auto& col : rho) {
      for (auto& row : col) {
        output_max = std::max(output_max, row);
      }
    }

    auto output_gain =
        static_cast<double>(std::numeric_limits<uint8_t>::max()) / output_max;

    for (int y = 0; y < n_pixels; ++y) {
      uint8_t row[n_pixels];
      for (auto x = 0; x < n_pixels; ++x) {
        row[x] = std::numeric_limits<uint8_t>::max() - output_gain * rho[y][x];
      }
      fw.write_row(row);
    }
  }

  template <class FileWriter>
  void output_kernel_bitmap(FileWriter& fw, F y, F angle) {

    F sec, sec_sq, A, B, C, bb_y, bb_z;
    detector.ellipse_bb(angle, 0, sec, sec_sq, A, B, C, bb_y, bb_z);

    Point emission_point = Point(0, 0);

    Pixel center_pixel =
        pixel_location(emission_point.first, emission_point.second);

    Pixel ur = Pixel(center_pixel.first - n_pixels_in_line(bb_y),
                     center_pixel.second + n_pixels_in_line(bb_z));
    Pixel dl = Pixel(center_pixel.first + n_pixels_in_line(bb_y),
                     center_pixel.second - n_pixels_in_line(bb_z));

    std::vector<std::pair<Pixel, F>> ellipse_kernels;
    ellipse_kernels.reserve(2000);

    std::vector<std::vector<F>> kernel_space;
    kernel_space.assign(bb_y, std::vector<F>(bb_z, 0));

    for (int iz = dl.second; iz < ur.second; ++iz) {
      for (int iy = ur.first; iy < dl.first; ++iy) {

        Point point = detector.pixel_center(iy, iz);

        if (detector.in_ellipse(A, B, C, emission_point, point)) {
          point -= emission_point;

          F event_kernel =
              kernel.calculate_kernel(
                  y, tan, sec, sec_sq, point, detector, sqrt_det_cor_mat) /
              sensitivity[iy][iz];
        }
      }
    }

    fw.template write_header<>(bb_y, bb_z);

    F output_max = 0;
    for (auto& col : kernel_space) {
      for (auto& row : col) {
        output_max = std::max(output_max, row);
      }
    }

    auto output_gain =
        static_cast<double>(std::numeric_limits<uint8_t>::max()) / output_max;

    for (int y = 0; y < n_pixels; ++y) {
      uint8_t row[n_pixels];
      for (auto x = 0; x < n_pixels; ++x) {
        row[x] = std::numeric_limits<uint8_t>::max() -
                 output_gain * kernel_space[y][x];
      }
      fw.write_row(row);
    }
  }

 private:
  int n_pixels_in_line(F length) const {
    return static_cast<int>((length + F(0.5)) / pixel_size);
  }

  void bb_update(Point ellipse_center, F angle, F y, F tan, int tid) {

    F sec, sec_sq, A, B, C, bb_y, bb_z;
    detector.ellipse_bb(angle, tan, sec, sec_sq, A, B, C, bb_y, bb_z);

    Pixel center_pixel =
        detector.pixel_location(ellipse_center.x, ellipse_center.y);

    Pixel ur = Pixel(center_pixel.x - n_pixels_in_line(bb_y),
                     center_pixel.y + n_pixels_in_line(bb_z));
    Pixel dl = Pixel(center_pixel.x + n_pixels_in_line(bb_y),
                     center_pixel.y - n_pixels_in_line(bb_z));

    std::vector<std::pair<Pixel, F>> ellipse_kernels;
    ellipse_kernels.reserve(2000);

    F acc = 0;
    for (int iz = dl.y; iz < ur.y; ++iz) {
      for (int iy = ur.x; iy < dl.x; ++iy) {
        Point point = detector.pixel_center(iy, iz);

        if (detector.in_ellipse(A, B, C, ellipse_center, point)) {
          point -= ellipse_center;

          F event_kernel = kernel(y,
                                  tan,
                                  sec,
                                  sec_sq,
                                  detector.radius,
                                  point,
                                  detector.inv_cor_mat_diag,
                                  sqrt_det_cor_mat) /
                           sensitivity[iy][iz];

          ellipse_kernels.push_back(
              std::make_pair(Pixel(iy, iz), event_kernel));
          acc += event_kernel * sensitivity[iy][iz] * rho[iy][iz];
        }
      }
    }

    for (auto& e : ellipse_kernels) {
      thread_rho[tid][e.first.x + (e.first.y * n_pixels)] +=
          e.second * rho[e.first.x][e.first.y] / acc;
    }
  }

  void simple_update(Point ellipse_center, F y, F z, int tid) {

    Pixel pixel = detector.pixel_location(ellipse_center);

    int y_line = 3 * (detector.sigma_z / pixel_size);
    int z_line = 3 * (detector.sigma_dl / pixel_size);

    std::vector<std::pair<Pixel, F>> ellipse_kernels;
    ellipse_kernels.reserve(2000);

    F acc = 0;

    for (int iz = pixel.x - y_line; iz < pixel.y + y_line; ++iz) {
      for (int iy = pixel.y - z_line; iy < pixel.z + z_line; ++iy) {

        Point point = detector.pixel_center(iy, iz);

        F event_kernel = kernel.test_kernel(
            y, z, point, detector.sigma_z, detector.sigma_dl);

        ellipse_kernels.push_back(std::make_pair(Pixel(iy, iz), event_kernel));
        acc += event_kernel * rho[iy][iz];
      }
    }

    for (auto& e : ellipse_kernels) {
      thread_rho[tid][e.first.first + (e.first.second * n_pixels)] +=
          e.second * rho[e.first.first][e.first.second] / acc;
    }
  }
};
