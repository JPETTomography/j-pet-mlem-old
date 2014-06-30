#pragma once

#include <cmath>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <unordered_map>
#include "event.h"

#include "util/png_writer.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"

#include "strip_detector.h"
#include "kernel.h"

#if _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#endif

#define RECONSTRUCTION_OLD_KERNEL 1

template <typename FType = double, typename D = StripDetector<FType>>
class Reconstruction {
 public:
  typedef FType F;
  typedef typename D::Pixel Pixel;
  typedef typename D::Point Point;
  F sqrt_det_correlation_matrix;

 private:
  static constexpr const F INVERSE_PI = Kernel<F>::INVERSE_PI;
  static constexpr const F INVERSE_POW_TWO_PI = Kernel<F>::INVERSE_POW_TWO_PI;

  int iteration;
  int n_pixels;
  F pixel_size;
  F pow_sigma_z, pow_sigma_dl;
  F inv_pow_sigma_z;
  F inv_pow_sigma_dl;

  std::vector<Event<F>> event_list;
  std::vector<std::vector<F>> rho;
  std::vector<std::vector<F>> rho_temp;
  std::vector<F> acc_log;
  std::vector<std::vector<F>> thread_rho;
  std::vector<std::vector<F>> lookup_table;

  D detector_;
  Kernel<F> kernel_;
  std::string output_;

 public:
  Reconstruction(int iteration, const D& detector)
      : iteration(iteration), detector_(detector) {
    init();
  };
  Reconstruction(int iteration,
                 F R_distance_a,
                 F scintilator_length,
                 int n_pixels,
                 F pixel_size,
                 F sigma_z,
                 F sigma_dl)
      : iteration(iteration),
        n_pixels(n_pixels),
        pixel_size(pixel_size),
        detector_(R_distance_a,
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
    kernel_ = Kernel<F>();
    rho.assign(n_pixels, std::vector<F>(n_pixels, F(100)));
    rho_temp.assign(n_pixels, std::vector<F>(n_pixels, F(10)));
    pow_sigma_z = detector_.sigma_z * detector_.sigma_z;
    pow_sigma_dl = detector_.sigma_dl * detector_.sigma_dl;

    sqrt_det_correlation_matrix = std::sqrt(
        detector_.inv_c(0, 0) * detector_.inv_c(1, 1) * detector_.inv_c(2, 2));

    inv_pow_sigma_z = F(1.0) / pow_sigma_z;
    inv_pow_sigma_dl = F(1.0) / pow_sigma_dl;

    lookup_table.assign(n_pixels, std::vector<F>(n_pixels));

    for (int y = 0; y < n_pixels; ++y) {
      for (int z = 0; z < n_pixels; ++z) {
        Point pp = pixel_center(y, z);
        lookup_table[y][z] = detector_.sensitivity(pp.first, pp.second);
      }
    }
  }

 public:
  /// Performs n_iterations of the list mode MEML algorithm
  void iterate(int n_iterations) {

    for (int i = 0; i < n_iterations; i++) {

      thread_rho.assign(omp_get_max_threads(),
                        std::vector<F>(n_pixels * n_pixels, F(0.0f)));

      std::cout << "ITERATION: " << i << std::endl;

      int size = event_list.size();

#if _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
      for (int id = 0; id < size; ++id) {

        int tid = omp_get_thread_num();

#if RECONSTRUCTION_OLD_KERNEL

        auto event = event_list[id];
        F tan = event.tan(detector_.radius);
        F y = event.y(tan);
        F z = event.z(y, tan);

        F angle = std::atan(tan);

        Point ellipse_center = Point(y, z);

        bb_pixel_updates(ellipse_center, angle, y, tan, tid);

#else

        F y = event_list[id].z_u;
        F z = event_list[id].z_d;

        if (id == 0) {

          std::cout << y << " " << z << std::endl;
        }

        Point ellipse_center = Point(y, z);

        simple_update(ellipse_center, y, z, tid, id);

#endif
      }

      rho.assign(n_pixels, std::vector<F>(n_pixels, F(0)));

      for (int i = 0; i < n_pixels; ++i) {
        for (int j = 0; j < n_pixels; ++j) {
          for (int k = 0; k < omp_get_max_threads(); ++k) {

            rho[i][j] += thread_rho[k][i + j * n_pixels];

            if (rho[i][j] > 0) {
              // std::cout << i << " " << j << " " << rho[i][j] << std::endl;
            }
          }
        }
      }
    }
  }

  void operator()(int n_blocks) { operator()(n_blocks, 1); }

  void operator()(int n_blocks, int n_iterations_in_block) {

    // F tan, y, z, angle;
    // Point ellipse_center;
    // Pixel pp;
    // int tid;

    std::cout << "Number of threads: " << omp_get_max_threads() << std::endl;

    for (int i = 0; i < n_blocks; i++) {
      std::cout << "ITERATION BLOCK: " << i << std::endl;

      iterate(n_iterations_in_block);
      // output reconstruction PNG

      std::string file = output_;
      file.append("_");
      file.append(std::to_string(i + 1));
      file.append(".png");

      png_writer png(file);
      png.write_header<>(n_pixels, n_pixels);

      F output_max = 0.0;
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
          row[x] =
              std::numeric_limits<uint8_t>::max() - output_gain * rho[y][x];
        }
        png.write_row(row);
      }
    }

    std::ofstream file;
    file.open("pixels_output.txt");
    for (int x = 0; x < n_pixels; ++x) {
      for (int y = 0; y < n_pixels; ++y) {

        if (rho[x][y] == 100) {
          rho[x][y] = 1;
        }

        file << x << " " << y << " " << rho[x][y] << std::endl;
      }
    }
  }

  void simple_update(Point& ellipse_center, F& y, F& z, int& tid, int& iter) {

    Pixel center_pixel =
        pixel_location(ellipse_center.first, ellipse_center.second);

    int y_line = 3 * (detector_.s_z() / pixel_size);
    int z_line = 3 * (detector_.s_dl() / pixel_size);

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

    F acc = F(0.0);

    for (int iz = center_pixel.first - y_line; iz < center_pixel.first + y_line;
         ++iz) {
      for (int iy = center_pixel.second - z_line;
           iy < center_pixel.second + z_line;
           ++iy) {

        Point pp = pixel_center(iy, iz);

        F event_kernel =
            kernel_.test_kernel(y, z, pp, detector_.s_z(), detector_.s_dl());

        ellipse_kernels.push_back(
            std::pair<Pixel, F>(Pixel(iy, iz), event_kernel));

        if (iter == 0) {

          std::cout << "PP:: " << iz << " " << iy << " " << pp.first << " "
                    << pp.second << " EVENT: " << event_kernel << std::endl;

          std::cout << "LOCATION: " << iy + (iz * n_pixels) << std::endl;
        }

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

    F cos_ = std::cos((angle));

    F sec_ = F(1.0) / cos_;
    F sec_sq_ = sec_ * sec_;

    F A = (((F(4.0) / (cos_ * cos_)) * inv_pow_sigma_dl) +
           (F(2.0) * tg * tg * inv_pow_sigma_z));
    F B = -F(4.0) * tg * inv_pow_sigma_z;
    F C = F(2.0) * inv_pow_sigma_z;
    F B_2 = (B / F(2.0)) * (B / F(2.0));
#if 0
    printf("A:= %f B:= %f c:= %f B_2:= %f\n", A, B, C, B_2);
#endif
    F bb_y = bby(A, C, B_2);

    F bb_z = bbz(A, C, B_2);
#if 0
    printf("bb_y:= %f bb_z:= %f\n", bb_y, bb_z);
#endif
    Pixel center_pixel =
        pixel_location(ellipse_center.first, ellipse_center.second);
#if 0
    printf("Center_Pixel y:= %d z:= %d\n",
           center_pixel.first,
           center_pixel.second);

#endif
    Pixel ur = Pixel(center_pixel.first - pixels_in_line(bb_y),
                     center_pixel.second + pixels_in_line(bb_z));
    Pixel dl = Pixel(center_pixel.first + pixels_in_line(bb_y),
                     center_pixel.second - pixels_in_line(bb_z));

    std::vector<std::pair<Pixel, F>> ellipse_kernels;
    ellipse_kernels.reserve(2000);
#if 0
    printf("iz:= %d Limit:= %d \n", dl.second, ur.first);
    printf("iy:= %d Limit:= %d \n", ur.first, dl.first);
#endif
    F acc = F(0.0);
    for (int iz = dl.second; iz < ur.second; ++iz) {
      for (int iy = ur.first; iy < dl.first; ++iy) {

        Point pp = pixel_center(iy, iz);

        if (in_ellipse(A, B, C, ellipse_center, pp)) {
#if 0
          printf("Pixel(%d,%d): %f %f\n", iy, iz, pp.first, pp.second);
#endif
          pp.first -= ellipse_center.first;
          pp.second -= ellipse_center.second;
#if 0
          printf("Pixel(%d,%d): SUB: %f %f\n", iy, iz, pp.first, pp.second);
#endif
          F event_kernel =
              kernel_.calculate_kernel(y,
                                       tg,
                                       sec_,
                                       sec_sq_,
                                       pp,
                                       detector_,
                                       sqrt_det_correlation_matrix) /
              lookup_table[iy][iz];

          ellipse_kernels.push_back(
              std::pair<Pixel, F>(Pixel(iy, iz), event_kernel));
          acc += event_kernel * lookup_table[iy][iz] * rho[iy][iz];
        }
      }
    }

    for (auto& e : ellipse_kernels) {

      thread_rho[tid][e.first.first + (e.first.second * n_pixels)] +=
          e.second * rho[e.first.first][e.first.second] / acc;
    }
  }

  bool in_ellipse(F& A, F& B, F& C, Point ellipse_center, Point p) {

    F dy = (p.first - ellipse_center.first);
    F dz = (p.second - ellipse_center.second);

    return (((A * (dy * dy)) + (B * dy * dz) + (C * (dz * dz)))) <= F(9.0)
               ? true
               : false;
  }

  F bbz(F& A, F& C, F& B_2) { return F(3.0) / std::sqrt(C - (B_2 / A)); }

  F bby(F& A, F& C, F& B_2) { return F(3.0) / std::sqrt(A - (B_2 / C)); }

  // coord Plane
  Pixel pixel_location(F y, F z) { return detector_.pixel_location(y, z); }

  // pixel Plane
  Point pixel_center(F y, F z) { return detector_.pixel_center(y, z); }

  int pixels_in_line(F length) {
    F d = (length + 0.5) / pixel_size;
    return int(d);
  }

 public:
  template <typename StreamType> Reconstruction& operator<<(StreamType& in) {

    int size;
    in >> size;
#if ___DEBUG_OUTPUT_SAVE
    std::cout << number_of_pixels << " " << pixel_s << " " << iter << " "
              << number_of_event_in_file << std::endl;
    std::cout << "VECTOR SIZE: " << event_list.size() << std::endl;

#endif
    for (int it = 0; it < size; ++it) {
      F z_u, z_d, dl;
      in >> z_u >> z_d >> dl;
      Event<F> temp_event(z_u, z_d, dl);
      event_list.push_back(temp_event);
    }
    return *this;
  }

  std::vector<Event<F>>& get_event_list() { return event_list; }

// FIXME: this confuses ICC
#ifndef __ICC
// template <typename StreamType>
// friend StreamType& operator>>(StreamType& in, Reconstruction& r) {
//  r << in;
//  return in;
// }
#endif

  void set_output_name(std::string output) { output_ = output; }
};
