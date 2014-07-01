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

#define RECONSTRUCTION_OLD_KERNEL 1

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
  Reconstruction(const Detector& detector) : detector(detector) { init(); }

  Reconstruction(F R_distance,
                 F scintilator_length,
                 int n_pixels,
                 F pixel_size,
                 F sigma_z,
                 F sigma_dl)
      : n_pixels(n_pixels),
        pixel_size(pixel_size),
        detector(R_distance,
                 scintilator_length,
                 n_pixels,
                 n_pixels,
                 pixel_size,
                 pixel_size,
                 sigma_z,
                 sigma_dl) {
    init();
  }

 private:
  void init() {
    kernel = Kernel<F>();
    rho.assign(n_pixels, std::vector<F>(n_pixels, 100));
    rho_temp.assign(n_pixels, std::vector<F>(n_pixels, 10));

    sqrt_det_cor_mat =
        std::sqrt(detector.inv_cor_mat_diag[0] * detector.inv_cor_mat_diag[1] *
                  detector.inv_cor_mat_diag[2]);

    sensitivity.assign(n_pixels, std::vector<F>(n_pixels));

    for (int y = 0; y < n_pixels; ++y) {
      for (int z = 0; z < n_pixels; ++z) {
        Point point = pixel_center(y, z);
        sensitivity[y][z] = detector.sensitivity(point.x, point.y);
      }
    }
  }

 public:
  /// Performs n_iterations of the list mode MEML algorithm
  template <typename ProgressCallback>
  void operator()(ProgressCallback progress,
                  int n_iterations,
                  int n_iterations_so_far) {

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

#if RECONSTRUCTION_OLD_KERNEL
        auto event = events[id];
        F tan = event.tan(detector.radius);
        F y = event.y(tan);
        F z = event.z(y, tan);

        F angle = std::atan(tan);

        Point ellipse_center = Point(y, z);

        bb_pixel_updates(ellipse_center, angle, y, tan, tid);
#else
        F y = events[id].z_u;
        F z = events[id].z_d;

        if (id == 0) {
          std::cout << y << " " << z << std::endl;
        }

        Point ellipse_center = Point(y, z);

        simple_update(ellipse_center, y, z, tid, id);
#endif
      }

      rho.assign(n_pixels, std::vector<F>(n_pixels, 0));

      for (int i = 0; i < n_pixels; ++i) {
        for (int j = 0; j < n_pixels; ++j) {
          for (int k = 0; k < omp_get_max_threads(); ++k) {

            rho[i][j] += thread_rho[k][i + j * n_pixels];
#if DEBUG
            if (rho[i][j] > 0) {
              std::cout << i << " " << j << " " << rho[i][j] << std::endl;
            }
#endif
          }
        }
      }
    }
  }

  void kernel_image(F& angle) {

    F y = 500;
    F z = 500;

    F tg = 0;

    F cos = std::cos((angle));

    F sec = 1 / cos;
    F sec_sq = sec * sec;

    F A = (4 / (cos * cos)) * inv_pow_sigma_dl + 2 * tg * tg * inv_pow_sigma_z;
    F B = -4 * tg * inv_pow_sigma_z;
    F C = 2 * inv_pow_sigma_z;
    F B_2 = (B / 2) * (B / 2);

    F bb_y = bby(A, C, B_2);
    F bb_z = bbz(A, C, B_2);

    Point emision_center = Point(0, 0);

    Pixel center_pixel =
        pixel_location(emision_center.first, emision_center.second);

    Pixel ur = Pixel(center_pixel.first - pixels_in_line(bb_y),
                     center_pixel.second + pixels_in_line(bb_z));
    Pixel dl = Pixel(center_pixel.first + pixels_in_line(bb_y),
                     center_pixel.second - pixels_in_line(bb_z));

    std::vector<std::pair<Pixel, F>> ellipse_kernels;
    ellipse_kernels.reserve(2000);

    std::vector<std::vector<F>> kernel_space;
    kernel_space.assign(bb_y, std::vector<F>(bb_z, 0));

    for (int iz = dl.second; iz < ur.second; ++iz) {
      for (int iy = ur.first; iy < dl.first; ++iy) {

        Point pp = pixel_center(iy, iz);

        if (in_ellipse(A, B, C, emision_center, pp)) {

          pp.first -= emision_center.first;
          pp.second -= emision_center.second;

          F event_kernel =
              kernel.calculate_kernel(
                  y, tg, sec, sec_sq, pp, detector, sqrt_det_cor_mat) /
              sensitivity[iy][iz];
        }
      }
    }

    png_writer png("kernel_space.png");
    png.write_header<>(bb_y, bb_z);

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
      png.write_row(row);
    }
  }

  void simple_update(Point& ellipse_center, F& y, F& z, int& tid, int& iter) {

    Pixel center_pixel =
        pixel_location(ellipse_center.first, ellipse_center.second);

    int y_line = 3 * (detector.s_z() / pixel_size);
    int z_line = 3 * (detector.s_dl() / pixel_size);

    if (iter == 0) {

      std::cout << "Pixel Limit: " << center_pixel.first << " "
                << center_pixel.second << " " << y_line << " " << z_line
                << std::endl;
      std::cout << "Limit min Limit max: " << center_pixel.first - y_line << " "
                << center_pixel.first + y_line << " "
                << center_pixel.second - z_line << " "
                << center_pixel.second + z_line << std::endl;
    }

    std::vector<std::pair<Pixel, F>> ellipse_kernels;
    ellipse_kernels.reserve(2000);

    F acc = 0;

    for (int iz = center_pixel.first - y_line; iz < center_pixel.first + y_line;
         ++iz) {
      for (int iy = center_pixel.second - z_line;
           iy < center_pixel.second + z_line;
           ++iy) {

        Point pp = pixel_center(iy, iz);

        F event_kernel =
            kernel.test_kernel(y, z, pp, detector.s_z(), detector.s_dl());

        ellipse_kernels.push_back(
            std::pair<Pixel, F>(Pixel(iy, iz), event_kernel));
#if DEBUG
        if (iter == 0) {

          std::cout << "PP:: " << iz << " " << iy << " " << pp.first << " "
                    << pp.second << " EVENT: " << event_kernel << std::endl;

          std::cout << "LOCATION: " << iy + (iz * n_pixels) << std::endl;
        }
#endif
        acc += event_kernel * rho[iy][iz];
      }
    }

    for (auto& e : ellipse_kernels) {

      thread_rho[tid][e.first.first + (e.first.second * n_pixels)] +=
          e.second * rho[e.first.first][e.first.second] / acc;
    }
  }

  void bb_pixel_updates(Point& ellipse_center,
                        F& angle,
                        F& y,
                        F& tg,
                        int& tid) {

    // std::cout << "0 " << tg << " " << y << std::endl;

    F cos = std::cos((angle));

    F sec = 1 / cos;
    F sec_sq = sec * sec;

    F A = (4 / (cos * cos)) * detector.inv_pow_sigma_dl +
          2 * tg * tg * detector.inv_pow_sigma_z;
    F B = -4 * tg * detector.inv_pow_sigma_z;
    F C = 2 * detector.inv_pow_sigma_z;
    F B_2 = (B / 2) * (B / 2);

#if DEBUG
    printf("A:= %f B:= %f c:= %f B_2:= %f\n", A, B, C, B_2);
#endif

    F bb_y = bby(A, C, B_2);
    F bb_z = bbz(A, C, B_2);

#if DEBUG
    printf("bb_y:= %f bb_z:= %f\n", bb_y, bb_z);
#endif

    Pixel center_pixel = pixel_location(ellipse_center.x, ellipse_center.y);

#if DEBUG
    printf("Center_Pixel y:= %d z:= %d\n",
           center_pixel.first,
           center_pixel.second);
#endif

    Pixel ur = Pixel(center_pixel.x - pixels_in_line(bb_y),
                     center_pixel.y + pixels_in_line(bb_z));
    Pixel dl = Pixel(center_pixel.x + pixels_in_line(bb_y),
                     center_pixel.y - pixels_in_line(bb_z));

    std::vector<std::pair<Pixel, F>> ellipse_kernels;
    ellipse_kernels.reserve(2000);

#if DEBUG
    printf("iz:= %d Limit:= %d \n", dl.second, ur.first);
    printf("iy:= %d Limit:= %d \n", ur.first, dl.first);
#endif

    F acc = 0;
    for (int iz = dl.y; iz < ur.y; ++iz) {
      for (int iy = ur.x; iy < dl.x; ++iy) {

        Point point = pixel_center(iy, iz);

        if (in_ellipse(A, B, C, ellipse_center, point)) {
          point.x -= ellipse_center.x;
          point.y -= ellipse_center.y;

          F event_kernel = kernel(y,
                                  tg,
                                  sec,
                                  sec_sq,
                                  detector.radius,
                                  point,
                                  detector.inv_cor_mat_diag,
                                  sqrt_det_cor_mat) /
                           sensitivity[iy][iz];

          ellipse_kernels.push_back(
              std::pair<Pixel, F>(Pixel(iy, iz), event_kernel));
          acc += event_kernel * sensitivity[iy][iz] * rho[iy][iz];
        }
      }
    }

    for (auto& e : ellipse_kernels) {
      thread_rho[tid][e.first.x + (e.first.y * n_pixels)] +=
          e.second * rho[e.first.x][e.first.y] / acc;
    }
  }

  bool in_ellipse(F& A, F& B, F& C, Point ellipse_center, Point p) {

    F dy = p.x - ellipse_center.x;
    F dz = p.y - ellipse_center.y;

    return (((A * (dy * dy)) + (B * dy * dz) + (C * (dz * dz)))) <= 9;
  }

  F bbz(F& A, F& C, F& B_2) { return 3 / std::sqrt(C - (B_2 / A)); }

  F bby(F& A, F& C, F& B_2) { return 3 / std::sqrt(A - (B_2 / C)); }

  // coord Plane
  Pixel pixel_location(F y, F z) { return detector.pixel_location(y, z); }

  // pixel Plane
  Point pixel_center(F y, F z) { return detector.pixel_center(y, z); }

  int pixels_in_line(F length) {
    F d = (length + F(0.5)) / pixel_size;
    return int(d);
  }

 public:
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

  std::vector<Event<F>>& get_event_list() { return events; }
};
