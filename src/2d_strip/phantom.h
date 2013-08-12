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

  typedef std::pair<int, int> Pixel;
  typedef std::pair<T, T> Point;

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
  std::vector<event<T> > event_list;
  std::vector<scintillator<> > scientilator_list;
  std::vector<T> output;
  static constexpr const T PI_2 = 1.5707963;
  static constexpr const T radian = T(M_PI / 180);

public:
  Phantom(int &iteration, int n_pixels, T &pixel_size, T &R_distance,
          T &Scentilator_length, T &x, T &y, T &a, T &b, T &phi)
      : iteration(iteration), n_pixels(n_pixels), pixel_size(pixel_size),
        R_distance(R_distance), Scentilator_length(Scentilator_length), _a(a),
        _b(b), _x(x), _y(y), _phi(phi) {
    _sin = sin(T(_phi * radian));
    _cos = cos(T(_phi * radian));
    _inv_a2 = T(1) / (_a * _a);
    _inv_b2 = T(1) / (_b * _b);
    output.assign(n_pixels * n_pixels, T(0.0));
  }

  bool in(T x, T y) const {
    /*
        T dx = x - _x;
        T dy = y - _y;

        T r_x = dx * _cos + dy * _sin;
        T r_y = -dx * _sin + dy * _cos;

        T r2 = r_x * r_x * _inv_a2 + r_y * r_y * _inv_b2;
      */

    T dy = (y - _y);
    T dz = (x - _x);
    T d1 = (_cos * dz + _sin * dy);
    T d2 = (-_cos * dy + _sin * dz);

    // std::cout << "ELLIPSE VALUE: " << (d1 * d1 / pow_sigma_z) +
    //                                       (d2 * d2 / pow_sigma_dl) <<
    // std::endl;

    return ((d1 * d1 / (_a * _a)) + (d2 * d2 / (_b * _b))) <= T(1) ? true
                                                                   : false;

    // return r2 < 1.0;
  }

  void emit_event(int n_threads) {

    std::vector<std::vector<event<T> > > event_list_per_thread;
    event<T> temp_event;

    T ry, rz, rangle;
    T z_u, z_d, dl;

#if OMP
    omp_set_num_threads(n_threads);
#endif

    event_list_per_thread.resize(omp_get_max_threads());

    std::uniform_real_distribution<T> uniform_dist(0, 1);
    std::uniform_real_distribution<T> uniform_y(_y - _b, _y + _b);
    std::uniform_real_distribution<T> uniform_z(_x - _a, _x + _a);
    std::normal_distribution<T> normal_dist_dz(0, 10);
    std::normal_distribution<T> normal_dist_dl(0, 63);
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
      rangle = M_PI_2 * uniform_dist(rng_list[omp_get_thread_num()]);

      if (in(rz, ry)) {

        z_u = rz + (R_distance - ry) * tan(rangle) +
              normal_dist_dz(rng_list[omp_get_thread_num()]);
        z_d = rz - (R_distance + ry) * tan(rangle) +
              normal_dist_dz(rng_list[omp_get_thread_num()]);
        dl = -T(2) * sqrt(T(1) + (tan(rangle) * tan(rangle))) +
             normal_dist_dl(rng_list[omp_get_thread_num()]);

        /*
                z_u = rz + (R_distance - ry) * tan(rangle);
                z_d = rz - (R_distance + ry) * tan(rangle);
                dl = -T(2) * ry * sqrt(T(1) + (tan(rangle) * tan(rangle)));

        */
        if (std::abs(z_u) < (Scentilator_length / T(2)) &&
            std::abs(z_d) < (Scentilator_length / T(2))) {
          // std::cout << z_u << " " << z_d << " " << dl << std::endl;

          T tan = event_tan(z_u, z_d);
          T _y = event_y(dl, tan);
          T _z = event_z(z_u, z_d, _y, tan);

          Pixel p = pixel_location(_y, _z);

          output[mem_location(p.first, p.second)] += T(1);

          temp_event.z_u = z_u;
          temp_event.z_d = z_d;
          temp_event.dl = dl;
          //  std::cout << "POINT: " << rz << " " << ry << " " <<
          // std::tan(rangle) <<
          //               " AFTERGAUSS: " << _z << " " << _y << " " << tan <<
          // std::endl;
          event_list_per_thread[omp_get_thread_num()].push_back(temp_event);
        }
      }
    }

    for (signed i = 0; i < omp_get_max_threads(); ++i) {

      event_list.insert(event_list.end(), event_list_per_thread[i].begin(),
                        event_list_per_thread[i].end());
    }

    // output reconstruction PNG
    png_writer png("phantom.png");
    png.write_header<>(n_pixels, n_pixels);

    double output_max = 0.0;
    for (auto it = output.begin(); it != output.end(); ++it) {
      output_max = std::max(output_max, *it);
    }

    auto output_gain =
        static_cast<double>(std::numeric_limits<uint8_t>::max()) / output_max;

    for (int y = 0; y < n_pixels; ++y) {
      uint8_t row[n_pixels];
      for (auto x = 0; x < n_pixels; ++x) {
        row[x] = std::numeric_limits<uint8_t>::max() -
                 output_gain * output[y * n_pixels + x];
      }
      png.write_row(row);
    }
  }

  T event_tan(T &z_u, T &z_d) const {
    return (z_u - z_d) / (T(2) * R_distance);
  }
  T event_y(T &dl, T &tan_event) const {
    return -T(0.5) * (dl / sqrt(T(1) + (tan_event * tan_event)));
  }
  T event_z(T &z_u, T &z_d, T &y, T &tan_event) const {
    return T(0.5) * (z_u + z_d + (T(2.0) * y * tan_event));
  }

  // coord Plane
  Pixel pixel_location(T y, T z) {
    return Pixel(std::floor((R_distance + y) / pixel_size),
                 std::floor((R_distance + z) / pixel_size));
  }

  int mem_location(int y, int z) { return int(y * n_pixels + z); }

  // pixel Plane
  Point pixel_center(T y, T z) {
    return Point(
        (std::floor((y) * pixel_size - R_distance)) + (T(0.5) * pixel_size),
        (std::floor((z) * pixel_size - R_distance)) + (T(0.5) * pixel_size));
  }

  template <typename StreamType> Phantom &operator>>(StreamType &out) {

    unsigned int n_pix = n_pixels;
    float pixel_s = pixel_size;
    unsigned int iter = iteration;
    unsigned int size = event_list.size();
    std::cout << "SAVE: " << n_pixels << " " << event_list.size() << std::endl;
    out << n_pix;
    out << pixel_s;
    out << iter;
    out << size;

    typename std::vector<event<T> >::iterator it;
    int i = 0;
    for (it = event_list.begin(); it != event_list.end(); ++it) {
      ++i;
      out << it->z_u << it->z_d << it->dl;
    }
    return *this;
  }

  template <typename StreamType>
  friend StreamType &operator<<(StreamType &out, Phantom &p) {
    p >> out;
    return out;
  }
};

#endif // STRIP_PET_H
