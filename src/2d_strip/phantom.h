#ifndef STRIP_PET_H
#define STRIP_PET_H

#include <cmath>
#include <vector>
#include <fstream>
#include <ctime>
#include <random>

#include <omp.h>

#include "data_structures.h"
#include "scintillator.h"
#include "bstream.h"

#define PI_2 1.5707963

typedef std::minstd_rand0 rng;

template <typename T = float>
class phantom {

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
  std::vector<event<T> > event_list;
  std::vector<scintillator<> > scientilator_list;

 public:
  phantom(int &iteration, int n_pixels, T &pixel_size, T &R_distance,
          T &Scentilator_length, T &x, T &y, T &a, T &b, T &phi)
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
    _sin = sin(static_cast<double>(_phi));
    _cos = cos(static_cast<double>(_phi));
    _inv_a2 = 1.0 / (_a * _a);
    _inv_b2 = 1.0 / (_b * _b);

    scientilator_list.push_back(
        scintillator<>(R_distance, 0.0f, Scentilator_length));
    scientilator_list.push_back(
        scintillator<>(-R_distance, 0.0f, Scentilator_length));
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

    std::vector<std::vector<event<T> > > event_list_per_thread;
    event<T> temp_event;

    T ry, rz, rangle;
    T z_u, z_d, dl;

    omp_set_num_threads(n_threads);

    event_list_per_thread.resize(omp_get_max_threads());

    std::uniform_real_distribution<T> uniform_dist(0, 1);
    std::normal_distribution<T> normal_dist(0, 10);

    rng rd;

    std::vector<rng> rng_list;

    for (int i = 0; i < omp_get_max_threads(); ++i) {

      rng_list.push_back(rd);
      rng_list[i].seed(42 + (3453 * i));
      // OR
      // Turn on leapfrogging with an offset that depends on the task id
    }

#pragma omp for schedule(static) private(ry, rz, rangle, z_u, z_d, dl)
    for (int i = 0; i < iteration; ++i) {

      ry = uniform_dist(rng_list[omp_get_thread_num()]);
      rz = uniform_dist(rng_list[omp_get_thread_num()]);
      rangle = PI_2 * uniform_dist(rng_list[omp_get_thread_num()]);

      if (in(ry, rz)) {

        z_u = rz + (R_distance - ry) * tan(rangle) +
              normal_dist(rng_list[omp_get_thread_num()]);
        z_d = rz - (R_distance + ry) * tan(rangle) +
              normal_dist(rng_list[omp_get_thread_num()]);
        dl = -static_cast<T>(2) *
                 sqrt(static_cast<T>(1) + (tan(rangle) * tan(rangle))) +
             normal_dist(rng_list[omp_get_thread_num()]);

        if (std::abs(z_u) < Scentilator_length &&
            std::abs(z_d) < Scentilator_length) {

          temp_event.z_u = z_u;
          temp_event.z_d = z_d;
          temp_event.dl = dl;

          event_list_per_thread[omp_get_thread_num()].push_back(temp_event);
        }
      }
    }

    for (signed i = 0; i < omp_get_max_threads(); ++i) {

      event_list.insert(event_list.end(), event_list_per_thread[i].begin(),
                        event_list_per_thread[i].end());
    }

    std::cout << "WIELKOSC: " << event_list.size() << std::endl;
  }

  pixel_location in_pixel(T &y, T &z) {

    return std::make_pair(((Scentilator_length + y) / pixel_size),
                          ((Scentilator_length + z) / pixel_size));
  }

  /*
  FileInt in_magic;
  in >> in_magic;
  if (in_magic != MAGIC_VERSION_TRIANGULAR &&
      in_magic != MAGIC_VERSION_FULL &&
      in_magic != MAGIC_VERSION_TOF_TRIANGULAR &&
      in_magic != MAGIC_VERSION_TOF_FULL && in_magic != MAGIC_VERSION_1 &&
      in_magic != MAGIC_VERSION_2) {
    throw("invalid file type format");
  }

  bool in_is_triangular =
      (in_magic != MAGIC_VERSION_FULL && in_magic != MAGIC_VERSION_TOF_FULL);
  bool in_is_tof = (in_magic == MAGIC_VERSION_TOF_TRIANGULAR ||
                    in_magic == MAGIC_VERSION_TOF_FULL);

  FileInt in_n_pixels_in_row;
  in >> in_n_pixels_in_row;
  if (in_is_triangular)
    in_n_pixels_in_row *= 2;

  FileInt in_n_emissions = 0;
  in >> in_n_emissions;

  FileInt in_n_detectors = 0;
  in >> in_n_detectors;

  FileInt in_n_tof_positions = 1;
  if (in_is_tof) {
    in >> in_n_tof_positions;
  }
#if DEBUG
  std::cerr << "in_n_pixels_in_row " << in_n_pixels_in_row << std::endl;
  std::cerr << "in_n_emissions " << in_n_emissions << std::endl;
  std::cerr << "in_n_detectors " << in_n_detectors << std::endl;
  std::cerr << "in_n_tof_positions " << in_n_tof_positions << std::endl;
#endif

  triangular_ = in_is_triangular;
  n_tof_positions_ = in_n_tof_positions;
  n_pixels_in_row_ = in_n_pixels_in_row;
  n_pixels_in_row_half_ = in_n_pixels_in_row / 2;
  n_emissions_ = in_n_emissions;
  n_detectors_ = in_n_detectors;
  n_2_detectors_ = in_n_detectors * 2;
  S n_detectors_2 = in_n_detectors / 2;
  S n_detectors_4 = in_n_detectors / 4;
  n_1_detectors_2_ = in_n_detectors + n_detectors_2;
  n_1_detectors_4_ = in_n_detectors + n_detectors_4;
  n_lors_ = LOR::end_for_detectors(n_detectors_).index();

  // load hits
  for (;;) {
    FileHalf a, b;
    in >> a >> b;
    LOR lor(a, b);

    if (in.eof())
      break;

    FileInt count;
    in >> count;

    // increment hits
    for (FileInt i = 0; i < count; ++i) {
      FileHalf x, y;
      FileInt position;
      FileInt hits;
      if (in_is_tof) {
        in >> position >> x >> y >> hits;
      } else {
        in >> x >> y >> hits;
        position = 0;
      }

      this->push_back(Element(lor, position, Pixel(x, y), hits));
    }
  }
}
*/

  void save_output(std::string fn) {

    obstream out(fn, std::ios::binary | std::ios::trunc);

    out << n_pixels;
    out << pixel_size;
    out << iteration;
    out << event_list.size();

    typename std::vector<event<T> >::iterator it;

    for (it = event_list.begin(); it != event_list.end(); ++it) {

      out << it->z_u << it->z_d << it->dl;
    }
  }

  void load_input(std::string fn) {

    ibstream in(fn, std::ios::binary);

    int number_of_pixels;
    float pixel_s;
    int iter;
    int number_of_event_in_file;

    in >> number_of_pixels;
    in >> pixel_s;
    in >> iter;
    in >> number_of_event_in_file;

    std::cout << number_of_pixels << " " << pixel_s << " " << iter << " "
              << number_of_event_in_file << std::endl;
  }
};

#endif  // STRIP_PET_H
