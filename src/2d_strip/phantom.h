#ifndef STRIP_PET_H
#define STRIP_PET_H

#include <cmath>
#include <vector>
#include <fstream>
#include <ctime>
#include <random>

#if OMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#endif

#include "event.h"
#include "scintillator.h"
#include "util/bstream.h"

typedef std::minstd_rand0 rng;

template <typename T = float> class Phantom {

  typedef std::pair<int, int> pixel_location;

 private:
  int iteration;
  int n_pixels;
  T pixel_size;
  T R_distance;
  T Scentilator_length;
  T _a, _b, _x, _y, _phi;
  T _sin;
  T _cos;
  T _inv_a2;
  T _inv_b2;
  std::vector<event<T>> event_list;
  std::vector<scintillator<>> scientilator_list;
  static constexpr const T PI_2 = 1.5707963;

 public:
  Phantom(int& iteration,
          int n_pixels,
          T& pixel_size,
          T& R_distance,
          T& Scentilator_length,
          T& x,
          T& y,
          T& a,
          T& b,
          T& phi)
      : iteration(iteration),
        n_pixels(n_pixels),
        pixel_size(pixel_size),
        R_distance(R_distance),
        Scentilator_length(Scentilator_length),
        _a(a),
        _b(b),
        _x(x),
        _y(y),
        _phi(phi) {
    _sin = sin(T(_phi));
    _cos = cos(T(_phi));
    _inv_a2 = T(1) / (_a * _a);
    _inv_b2 = T(1) / (_b * _b);
  }

  bool in(T x, T y) const {

    T dx = x - _x;
    T dy = y - _y;

    T r_x = dx * _cos + dy * _sin;
    T r_y = -dx * _sin + dy * _cos;

    T r2 = r_x * r_x * _inv_a2 + r_y * r_y * _inv_b2;

    return r2 < 1.0;
  }

  void emit_event(int n_threads) {

    std::vector<std::vector<event<T>>> event_list_per_thread;
    event<T> temp_event;

    T ry, rz, rangle;
    T z_u, z_d, dl;

#if OMP
    omp_set_num_threads(n_threads);
#endif

    event_list_per_thread.resize(omp_get_max_threads());

    std::uniform_real_distribution<T> uniform_dist(0, 1);
    std::uniform_real_distribution<T> uniform_y(_y - _a, _y + _a);
    std::uniform_real_distribution<T> uniform_z(_x - _b, _x + _b);
    std::normal_distribution<T> normal_dist(0, 10);
    rng rd;

    std::vector<rng> rng_list;

    for (int i = 0; i < omp_get_max_threads(); ++i) {

      rng_list.push_back(rd);
      rng_list[i].seed(42 + (3453 * i));
      // OR
      // Turn on leapfrogging with an offset that depends on the task id
    }
#if OMP
#pragma omp for schedule(static) private(ry, rz, rangle, z_u, z_d, dl)
#endif
    for (int i = 0; i < iteration; ++i) {

      ry = uniform_y(rng_list[omp_get_thread_num()]);
      rz = uniform_z(rng_list[omp_get_thread_num()]);
      rangle = PI_2 * uniform_dist(rng_list[omp_get_thread_num()]);

      if (in(ry, rz)) {

        z_u = rz + (R_distance - ry) * tan(rangle) +
              normal_dist(rng_list[omp_get_thread_num()]);
        z_d = rz - (R_distance + ry) * tan(rangle) +
              normal_dist(rng_list[omp_get_thread_num()]);
        dl = -T(2) * sqrt(T(1) + (tan(rangle) * tan(rangle))) +
             normal_dist(rng_list[omp_get_thread_num()]);

        if (std::abs(z_u) < (Scentilator_length / T(2)) &&
            std::abs(z_d) < (Scentilator_length / T(2))) {
          // std::cout << z_u << " " << z_d << " " << dl << std::endl;
          temp_event.z_u = z_u;
          temp_event.z_d = z_d;
          temp_event.dl = dl;

          event_list_per_thread[omp_get_thread_num()].push_back(temp_event);
        }
      }
    }

    for (signed i = 0; i < omp_get_max_threads(); ++i) {

      event_list.insert(event_list.end(),
                        event_list_per_thread[i].begin(),
                        event_list_per_thread[i].end());
    }
  }

  void save_output(std::string fn) {

    obstream out(fn, std::ios::binary | std::ios::trunc);

    unsigned int n_pix = n_pixels;
    float pixel_s = pixel_size;
    unsigned int iter = iteration;
    unsigned int size = event_list.size();
    out << n_pix;
    out << pixel_s;
    out << iter;
    out << size;

    typename std::vector<event<T>>::iterator it;
    int i = 0;
    for (it = event_list.begin(); it != event_list.end(); ++it) {
      ++i;
      out << it->z_u << it->z_d << it->dl;
    }
  }
};

#endif  // STRIP_PET_H
