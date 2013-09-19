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

#include"strip_detector.h"

#if OMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#endif



template <typename T = float> class Reconstruction {
 public:
  typedef std::pair<int, int> Pixel;
  typedef std::pair<T, T> Point;

 private:

  static constexpr const T INVERSE_PI = T(1.0 / M_PI);
  static constexpr const T INVERSE_POW_TWO_PI = T(1.0 / (2.0 * M_PI * M_PI));

  int iteration;
  T R_distance;
  T Scentilator_length;
  int n_pixels;
  T pixel_size;
  T sigma_z, sigma_dl;
  T pow_sigma_z, pow_sigma_dl;
  T inv_pow_sigma_z;
  T inv_pow_sigma_dl;
  T half_scentilator_length;

  std::vector<event<T>> event_list;
  std::vector<std::vector<T>> inverse_correlation_matrix;
  T sqrt_det_correlation_matrix;
  std::vector<std::vector<T>> rho;
  std::vector<std::vector<T>> rho_temp;
  std::vector<T> acc_log;
  std::vector<std::vector<T>> thread_rho;
  std::vector<std::vector<T>> lookup_table;

 public:
  Reconstruction(int iteration,
                 T R_distance,
                 T Scentilator_length,
                 int n_pixels,
                 T pixel_size,
                 T sigma_z,
                 T sigma_dl)
      : iteration(iteration),
        R_distance(R_distance),
        Scentilator_length(Scentilator_length),
        n_pixels(n_pixels),
        pixel_size(pixel_size),
        sigma_z(sigma_z),
        sigma_dl(sigma_dl) {

    rho.assign(n_pixels, std::vector<T>(n_pixels, T(100)));
    rho_temp.assign(n_pixels, std::vector<T>(n_pixels, T(10)));
    pow_sigma_z = sigma_z * sigma_z;
    pow_sigma_dl = sigma_dl * sigma_dl;

    inverse_correlation_matrix.resize(3, std::vector<T>(3, T()));

    inverse_correlation_matrix[0][0] = T(1) / pow_sigma_z;
    inverse_correlation_matrix[0][1] = T(0.0f);
    inverse_correlation_matrix[0][2] = T(0.0f);

    inverse_correlation_matrix[1][0] = T(0.0f);
    inverse_correlation_matrix[1][1] = T(1) / pow_sigma_z;
    inverse_correlation_matrix[1][2] = T(0.0f);

    inverse_correlation_matrix[2][0] = T(0.0f);
    inverse_correlation_matrix[2][1] = T(0.0f);
    inverse_correlation_matrix[2][2] = T(1) / pow_sigma_dl;

    sqrt_det_correlation_matrix = std::sqrt(inverse_correlation_matrix[0][0] *
                                            inverse_correlation_matrix[1][1] *
                                            inverse_correlation_matrix[2][2]);

    inv_pow_sigma_z = T(1.0) / pow_sigma_z;
    inv_pow_sigma_dl = T(1.0) / pow_sigma_dl;
    half_scentilator_length = Scentilator_length * T(0.5f);

    lookup_table.assign(n_pixels, std::vector<T>(n_pixels));

    for (int y = 0; y < n_pixels; ++y) {
      for (int z = 0; z < n_pixels; ++z) {

        Point pp = pixel_center(y, z);
        lookup_table[y][z] = sensitivity(pp.first, pp.second);
      }
    }
  }

  T multiply_elements(T* vec_a, T* vec_b) {

    T output = T(0.0);

    output += vec_a[0] * inverse_correlation_matrix[0][0] * vec_b[0];
    output += vec_a[1] * inverse_correlation_matrix[1][1] * vec_b[1];
    output += vec_a[2] * inverse_correlation_matrix[2][2] * vec_b[2];

    return output;
  }

  float fexp(float& x) {
    volatile union {
      T f;
      unsigned int i;
    } cvt;

    /* exp(x) = 2^i * 2^f; i = floor (log2(e) * x), 0 <= f <= 1 */
    T t = x * 1.442695041;
    T fi = floorf(t);
    T f = t - fi;
    int i = (int)fi;
    cvt.f =
        (0.3371894346f * f + 0.657636276f) * f + 1.00172476f; /* compute 2^f */
    cvt.i += (i << 23);                                       /* scale by 2^i */
    return cvt.f;
  }

  double fexp(double& x) { return std::exp(x); }

  T kernel(T& y, T& _tan, T& inv_cos, T& pow_inv_cos, Point& pixel_center) {

    T vec_o[3];
    T vec_a[3];
    T vec_b[3];

    vec_o[0] = -(pixel_center.first + y - R_distance) * _tan * pow_inv_cos;
    vec_o[1] = -(pixel_center.first + y + R_distance) * _tan * pow_inv_cos;
    vec_o[2] =
        -(pixel_center.first + y) * inv_cos * (T(1) + T(2) * (_tan * _tan));

    vec_a[0] = -(pixel_center.first + y - R_distance) * pow_inv_cos;
    vec_a[1] = -(pixel_center.first + y + R_distance) * pow_inv_cos;
    vec_a[2] = -T(2) * (pixel_center.first + y) * (inv_cos * _tan);

    vec_b[0] = pixel_center.second - (pixel_center.first * _tan);
    vec_b[1] = pixel_center.second - (pixel_center.first * _tan);
    vec_b[2] = -T(2) * pixel_center.first * inv_cos;

    T a_ic_a = multiply_elements(vec_a, vec_a);
    T b_ic_a = multiply_elements(vec_b, vec_a);
    T b_ic_b = multiply_elements(vec_b, vec_b);
    T o_ic_b = multiply_elements(vec_o, vec_b);

    T norm = a_ic_a + (T(2.f) * o_ic_b);

    T element_before_exp =
        INVERSE_POW_TWO_PI * (sqrt_det_correlation_matrix / std::sqrt(norm));

    T exp_element = -T(0.5f) * (b_ic_b - ((b_ic_a * b_ic_a) / norm));

    T _exp = fexp(exp_element);

    return (element_before_exp * _exp);
  }

  T sensitivity(T y, T z) {

    T L_plus = (half_scentilator_length + z);
    T L_minus = (half_scentilator_length - z);
    T R_plus = R_distance + y;
    T R_minus = R_distance - y;
    return INVERSE_PI *
           (std::atan(std::min(L_minus / R_minus, L_plus / R_plus)) -
            std::atan(std::max(-L_plus / R_minus, -L_minus / R_plus)));
  }

  /** Performs n_iterations of the list mode MEML algorithm
   */
  void iterate(int n_iterations) {

    for (int i = 0; i < n_iterations; i++) {

      thread_rho.assign(omp_get_max_threads(),
                        std::vector<T>(n_pixels * n_pixels, T(0.0f)));

      std::cout << "ITERATION: " << i << std::endl;

      int size = event_list.size();

#if OMP
#pragma omp parallel for schedule(dynamic)
#endif
      for (int id = 0; id < size; ++id) {

        int tid = omp_get_thread_num();

        T tan = event_tan(event_list[id].z_u, event_list[id].z_d, R_distance);
        T y = event_y(event_list[id].dl, tan);
        T z = event_z(event_list[id].z_u, event_list[id].z_d, y, tan);

        T angle = std::atan(tan);

        Point ellipse_center = Point(y, z);

        Pixel pp = pixel_location(y, z);

        bb_pixel_updates(ellipse_center, angle, y, tan, tid);
      }

      rho.assign(n_pixels, std::vector<T>(n_pixels, T(0)));

      for (int i = 0; i < n_pixels; ++i) {
        for (int j = 0; j < n_pixels; ++j) {
          for (int k = 0; k < omp_get_max_threads(); ++k) {

            rho[i][j] += thread_rho[k][i + j * n_pixels];
          }
        }
      }
    }
  }

  void operator()(int n_blocks) { operator()(n_blocks, 1); }
  void operator()(int n_blocks, int n_iterations_in_block) {

    // T tan, y, z, angle;
    // Point ellipse_center;
    // Pixel pp;
    // int tid;

    std::cout << "Number of threads: " << omp_get_max_threads() << std::endl;

    for (int i = 0; i < n_blocks; i++) {
      std::cout << "ITERATION BLOCK: " << i << std::endl;

      iterate(n_iterations_in_block);
      // output reconstruction PNG

      std::string file = std::string("rec_iteration_");

      file.append(std::to_string(i + 1));
      file.append(".png");

      png_writer png(file);
      png.write_header<>(n_pixels, n_pixels);

      T output_max = 0.0;
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

  void bb_pixel_updates(Point& ellipse_center,
                        T& angle,
                        T& y,
                        T& tg,
                        int& tid) {

    T acc = T(0.0);

    T _cos = std::cos((angle));

    T inv_cos = T(1.0) / _cos;
    T pow_inv_cos = inv_cos * inv_cos;

    T A = (((T(4.0) / (_cos * _cos)) * inv_pow_sigma_dl) +
           (T(2.0) * tg * tg * inv_pow_sigma_z));
    T B = -T(4.0) * tg * inv_pow_sigma_z;
    T C = T(2.0) * inv_pow_sigma_z;
    T B_2 = (B / T(2.0)) * (B / T(2.0));

    T bb_z = T(3.0) / std::sqrt(C - (B_2 / A));

    T bb_y = T(3.0) / std::sqrt(A - (B_2 / C));

    Pixel center_pixel =
        pixel_location(ellipse_center.first, ellipse_center.second);

    Pixel ur = Pixel(center_pixel.first - pixels_in_line(bb_y),
                     center_pixel.second + pixels_in_line(bb_z));
    Pixel dl = Pixel(center_pixel.first + pixels_in_line(bb_y),
                     center_pixel.second - pixels_in_line(bb_z));

    std::vector<std::pair<Pixel, T>> ellipse_kernels;
    ellipse_kernels.reserve(2000);

    Point pp = pixel_center(ur.first, dl.second);

    int count = 0;

    for (int iz = dl.second; iz < ur.second; ++iz) {
      for (int iy = ur.first; iy < dl.first; ++iy) {

        Point pp = pixel_center(iy, iz);

        if (in_ellipse(A, B, C, ellipse_center, pp)) {

          pp.first -= ellipse_center.first;
          pp.second -= ellipse_center.second;

          T event_kernel =
              kernel(y, tg, inv_cos, pow_inv_cos, pp) / lookup_table[iy][iz];

          ellipse_kernels.push_back(
              std::pair<Pixel, T>(Pixel(iy, iz), event_kernel));
          acc += event_kernel * lookup_table[iy][iz] * rho[iy][iz];

          count++;
        }
      }
    }
    for (auto& e : ellipse_kernels) {

      thread_rho[tid][e.first.first + (e.first.second * n_pixels)] +=
          e.second * rho[e.first.first][e.first.second] / acc;
    }
  }

  bool in_ellipse(T& A, T& B, T& C, Point ellipse_center, Point p) {

    T dy = (p.first - ellipse_center.first);
    T dz = (p.second - ellipse_center.second);

    return (((A * (dy * dy)) + (B * dy * dz) + (C * (dz * dz)))) <= T(9.0)
               ? true
               : false;
  }

  // coord Plane
  Pixel pixel_location(T y, T z) {

    return Pixel(std::floor((R_distance - y) / pixel_size),
                 std::floor((R_distance + z) / pixel_size));
  }

  // pixel Plane
  Point pixel_center(T y, T z) {

    return Point((R_distance - (y + T(0.5)) * pixel_size),
                 (z + T(0.5)) * pixel_size - R_distance);
  }

  int pixels_in_line(T& length) {
    T d = length / pixel_size;
    return (d > 0) ? int(d) : int(d + 1);
  }

  template <typename StreamType> Reconstruction& operator<<(StreamType& in) {

    event<T> temp_event;

    int i = 0;
#if DEBUG_OUTPUT_SAVE
    std::cout << number_of_pixels << " " << pixel_s << " " << iter << " "
              << number_of_event_in_file << std::endl;
    std::cout << "VECTOR SIZE: " << event_list.size() << std::endl;

#endif
    for (;;) {

      if (in.eof()) {
        break;
      }
      i++;

      T z_u, z_d, dl;

      in >> z_u >> z_d >> dl;

      temp_event.z_u = z_u;
      temp_event.z_d = z_d;
      temp_event.dl = dl;

      T tan = event_tan(z_u, z_d, R_distance);
      T y = event_y(z_d, tan);
      T z = event_z(z_u, z_d, y, tan);

      Pixel p = pixel_location(y, z);

      event_list.push_back(temp_event);
    }

    return *this;
  }

  template <typename StreamType>
  friend StreamType& operator>>(StreamType& in, Reconstruction& r) {
    r << in;
    return in;
  }
};

