#pragma once

#include <cmath>
#include <vector>
#include <ctime>
#include <random>
#include <algorithm>

#include "event.h"

#if _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#endif

#define MAIN_PHANTOM 1

template <typename F> int sgn(F val) { return (0 < val) - (val < 0); }

template <typename FType = double> class Phantom {
  typedef FType F;
  typedef std::pair<int, int> Pixel;
  typedef std::pair<F, F> Point;
  typedef std::minstd_rand0 rng;

 private:
  int n_pixels;
  F pixel_size;
  F R_distance;
  F scintillator_length;
  F sin;
  F cos;
  F inv_a2;
  F inv_b2;
  F sigma_z;
  F sigma_dl;
  std::vector<EllipseParameters<F>> ellipse_list;
  std::vector<Event<F>> events;
  std::vector<std::vector<F>> output;
  std::vector<std::vector<F>> output_without_errors;

 public:
  Phantom(std::vector<EllipseParameters<F>>& el,
          int n_pixels,
          F pixel_size,
          F R_distance,
          F Scentilator_length,
          F sigma_z,
          F sigma_dl)
      : n_pixels(n_pixels),
        pixel_size(pixel_size),
        R_distance(R_distance),
        scintillator_length(Scentilator_length),
        sigma_z(sigma_z),
        sigma_dl(sigma_dl) {
    ellipse_list = el;
    output.assign(n_pixels, std::vector<F>(n_pixels, 0));
    output_without_errors.assign(n_pixels, std::vector<F>(n_pixels, 0));
  }

  bool in(F y, F z, EllipseParameters<F> el) const {

    F dy = (y - el.y);
    F dz = (z - el.x);
    F d1 = (sin * dy + cos * dz);  // y
    F d2 = (sin * dz - cos * dy);  // z

    return (d1 * d1 / (el.a * el.a)) + (d2 * d2 / (el.b * el.b)) <= 1;
  }

  void operator()() {
    const F RADIAN = F(M_PI / 180);

    std::vector<std::vector<Event<F>>> event_list_per_thread;

    event_list_per_thread.resize(omp_get_max_threads());

    for (auto& el : ellipse_list) {

      F max = std::max(el.a, el.b);

      std::cout << el.x << " " << el.y << " " << el.a << " " << el.b << " "
                << el.angle << std::endl;

      std::uniform_real_distribution<F> uniform_angle(-1, 1);
      std::uniform_real_distribution<F> uniform_y(el.y - max, el.y + max);
      std::uniform_real_distribution<F> uniform_z(el.x - max, el.x + max);
      std::normal_distribution<F> normal_dist_dz(0, sigma_z);
      std::normal_distribution<F> normal_dist_dl(0, sigma_dl);
      rng rd;

      std::vector<rng> rng_list;

      for (int i = 0; i < omp_get_max_threads(); ++i) {

        rng_list.push_back(rd);
        rng_list[i].seed(42 + (3453 * i));
        // OR
        // Turn on leapfrogging with an offset that depends on the task id
      }

      sin = std::sin(el.angle * RADIAN);
      cos = std::cos(el.angle * RADIAN);
      inv_a2 = 1 / (el.a * el.a);
      inv_b2 = 1 / (el.b * el.b);

      int n_emissions = el.n_emissions;  // FIXME: el.n_emission is float
#if _OPENMP
#pragma omp for schedule(static)
#endif
      for (int emission = 0; emission < n_emissions; ++emission) {

#if MAIN_PHANTOM
        F ry = uniform_y(rng_list[omp_get_thread_num()]);
        F rz = uniform_z(rng_list[omp_get_thread_num()]);
        F rangle = F(M_PI_4) * uniform_angle(rng_list[omp_get_thread_num()]);

        if (in(ry, rz, el) /* && std::abs(rangle) != M_PI_2 */) {
          F z_u, z_d, dl;

          z_u = rz + (R_distance - ry) * std::tan(rangle);
          z_d = rz - (R_distance + ry) * std::tan(rangle);
          dl = -2 * ry * std::sqrt(1 + (std::tan(rangle) * std::tan(rangle)));

          z_u += normal_dist_dz(rng_list[omp_get_thread_num()]);
          z_d += normal_dist_dz(rng_list[omp_get_thread_num()]);
          dl += normal_dist_dl(rng_list[omp_get_thread_num()]);

          if (std::abs(z_u) < scintillator_length / 2 &&
              std::abs(z_d) < scintillator_length / 2) {

            Event<F> event(z_u, z_d, dl);

            F tan, y, z;
            event.transform(R_distance, tan, y, z);

            Pixel p = pixel_location(y, z);
            Pixel pp = pixel_location(ry, rz);

            output[p.first][p.second]++;
            output_without_errors[pp.first][pp.second]++;

            event_list_per_thread[omp_get_thread_num()].push_back(event);
          }
        }
#else
        F ry = el.y + normal_dist_dl(rng_list[omp_get_thread_num()]);
        F rz = el.x + normal_dist_dz(rng_list[omp_get_thread_num()]);

        F y = ry;
        F z = rz;

        Pixel p = pixel_location(x, y);
        Pixel pp = pixel_location(0, 0);

        output[p.first][p.second]++;
        output_without_errors[pp.first][pp.second]++;

        Event<F> event(ry, rz, dl);
        event_list_per_thread[omp_get_thread_num()].push_back(event);
#endif
      }

      for (int i = 0; i < omp_get_max_threads(); ++i) {
        events.insert(events.end(),
                      event_list_per_thread[i].begin(),
                      event_list_per_thread[i].end());
      }

      std::cout << "VECTOR SIZE: " << events.size() << std::endl;

    }
  }

  template <class FileWriter>
  void output_bitmap(FileWriter& fw, bool wo_errors = false) {
    fw.template write_header<>(n_pixels, n_pixels);

    auto& target_output = wo_errors ? output_without_errors : output;
    F output_max = 0;
    for (auto& col : target_output) {
      for (auto& row : col)
        output_max = std::max(output_max, row);
    }

    auto output_gain =
        static_cast<double>(std::numeric_limits<uint8_t>::max()) / output_max;

    uint8_t* row = (uint8_t*)alloca(n_pixels);
    for (int y = 0; y < n_pixels; ++y) {
      for (auto x = 0; x < n_pixels; ++x) {
        row[x] = std::numeric_limits<uint8_t>::max() -
                 output_gain * target_output[y][x];
      }
      fw.write_row(row);
    }
  }

  // coord Plane
  Pixel pixel_location(F y, F z) {
    return Pixel((R_distance - y) / pixel_size, (R_distance + z) / pixel_size);
  }

  template <typename StreamType> Phantom& operator>>(StreamType& out) {

    int size = events.size();
    out << size;
    for (auto& event : events) {
      out << event.z_u << event.z_d << event.dl;
    }
    return *this;
  }

  template <typename StreamType>
  friend StreamType& operator<<(StreamType& out, Phantom& p) {
    p >> out;
    return out;
  }
};
