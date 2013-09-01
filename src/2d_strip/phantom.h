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

template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }

template <typename T = double> class Phantom {

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
  std::vector<std::vector<T> > output;
  static constexpr const T PI_2 = 1.5707963;
  static constexpr const T radian = T(M_PI / 180);

public:
  Phantom(int &iteration, int n_pixels, T &pixel_size, T &R_distance,
          T &Scentilator_length, T &x, T &y, T &a, T &b, T &phi)
      : iteration(iteration), n_pixels(n_pixels), pixel_size(pixel_size),
        R_distance(R_distance), Scentilator_length(Scentilator_length), _a(a),
        _b(b), _x(x), _y(y), _phi(phi) {
    _sin = std::sin(T(_phi * radian));
    _cos = std::cos(T(_phi * radian));
    _inv_a2 = T(1) / (_a * _a);
    _inv_b2 = T(1) / (_b * _b);
    output.assign(n_pixels, std::vector<T>(n_pixels, T(0)));
  }

  bool in(T y, T z) const {

    T pow_sigma_dl = 63 * 63;
    T pow_sigma_z = 10 * 10;

    T dy = (y - _y);
    T dz = (z - _x);
    T d1 = (_sin * dy + _cos * dz); // y
    T d2 = (_sin * dz - _cos * dy); // x

    /*
          T tg = (_sin / _cos);

          T A = (((T(4.0) * (T(1.0) / (_cos * _cos))) / pow_sigma_dl) +
                  (T(2.0) * tg * tg / pow_sigma_z));
          T B = -T(4.0) * tg / pow_sigma_z;
          T C = T(2.0) / pow_sigma_z;

            return ((A * (d2 * d2)) + (B * d1 * d2) + (C * (d1 * d1))) <= T(1) ?
       true:false;
    */

    return ((d1 * d1 / (_a * _a)) + (d2 * d2 / (_b * _b))) <= T(1) ? true
                                                                   : false;
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
    std::uniform_real_distribution<T> uniform_z(_x - _b, _x + _b);
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
      rangle = (M_PI_4 - radian) * uniform_dist(rng_list[omp_get_thread_num()]);

      if (in(ry, rz)) {

        z_u = rz + (R_distance - ry) * tan(rangle);
        z_d = rz - (R_distance + ry) * tan(rangle);
        dl = -T(2) * ry * sqrt(T(1) + (tan(rangle) * tan(rangle)));

        /*
                z_u = rz + (R_distance - ry) * tan(rangle) +
                      normal_dist_dz(rng_list[omp_get_thread_num()]);
                z_d = rz - (R_distance + ry) * tan(rangle) +
                      normal_dist_dz(rng_list[omp_get_thread_num()]);
                dl = -T(2) * ry * sqrt(T(1) + (tan(rangle) * tan(rangle))) +
                     normal_dist_dl(rng_list[omp_get_thread_num()]);
        */
        if (std::abs(z_u) < (Scentilator_length / T(2)) &&
            std::abs(z_d) < (Scentilator_length / T(2))) {

          Pixel p = pixel_location(ry, rz);

          //  std::cout << p.first << " " << p.second << " " << n_pixels <<
          // std::endl;

          output[p.first][p.second]++;

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

    T output_max = 0.0;
    for (auto &col : output) {
      for (auto &row : col)
        output_max = std::max(output_max, row);
    }

    auto output_gain =
        static_cast<double>(std::numeric_limits<uint8_t>::max()) / output_max;

    for (int y = 0; y < n_pixels; ++y) {
      uint8_t row[n_pixels];
      for (auto x = 0; x < n_pixels; ++x) {
        row[x] =
            std::numeric_limits<uint8_t>::max() - output_gain * output[y][x];
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
    return Pixel(std::floor((R_distance - y) / pixel_size),
                 std::floor((R_distance + z) / pixel_size));
  }

  Point pixel_center(T y, T z) {

    int sgn_y = sgn<T>(y);
    int sgn_z = sgn<T>(z);

    return Point((std::ceil(R_distance - (y) * pixel_size)) +
                     (sgn<T>(y) * T(0.5) * pixel_size),
                 (std::ceil((z) * pixel_size - R_distance)) +
                     (sgn<T>(z) * T(0.5) * pixel_size));
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
