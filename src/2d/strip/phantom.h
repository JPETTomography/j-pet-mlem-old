#pragma once

#include <cmath>
#include <vector>
#include <ctime>
#include <random>
#include <algorithm>

#include "event.h"
#include "2d/geometry/ellipse.h"

#if _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#endif

namespace PET2D {
namespace Strip {

template <typename F> int sgn(F val) { return (0 < val) - (val < 0); }

/// Virtual phantom region made of ellipse and intensity
template <typename FType> struct PhantomRegion {
  using F = FType;

  PhantomRegion(const Ellipse<F>& ellipse, F intensity)
      : shape(ellipse), intensity(intensity), weight(intensity * shape.area) {}

  bool contains(Point<F> p) const { return shape.constains(p); }

  const Ellipse<F> shape;
  const F intensity;
  const F weight;
};

/// Virtual phantom made of elliptical regions
template <typename D, typename FType, typename SType> class Phantom {
  using F = FType;
  using S = SType;
  using Pixel = PET2D::Pixel<S>;
  using rng = std::minstd_rand0;

 private:
  D detector;

  std::vector<PhantomRegion<F>> region_list;
  std::vector<F> CDF;
  std::vector<EllipsePointGenerator<F>> point_generators;

  std::vector<Event<F>> events;
  std::vector<std::vector<F>> output;
  std::vector<std::vector<F>> output_without_errors;

  std::uniform_real_distribution<F> uniform;
  std::uniform_real_distribution<F> uniform_angle;

 public:
  Phantom(const D& detector, const std::vector<PhantomRegion<F>>& el)
      : detector(detector),
        region_list(el),
        CDF(el.size(), 0),
        uniform_angle(-1, 1) {
    CDF[0] = region_list[0].weight;

    for (size_t i = 1; i < el.size(); i++) {
      CDF[i] = region_list[i].weight + CDF[i - 1];
    }
    F norm = CDF[el.size() - 1];
    for (size_t i = 0; i < el.size(); i++) {
      CDF[i] /= norm;
    }

    for (size_t i = 0; i < el.size(); ++i)
      point_generators.emplace_back(el[i].shape);

    output.assign(detector.n_y_pixels, std::vector<F>(detector.n_z_pixels, 0));
    output_without_errors.assign(detector.n_y_pixels,
                                 std::vector<F>(detector.n_z_pixels, 0));
  }

  size_t n_events() { return events.size(); }

  template <typename G> size_t choose_region(G& gen) {
    F r = uniform(gen);
    size_t i = 0;

    while (r > CDF[i])
      ++i;

    return i;
  }

  template <typename Generator> Point<F> gen_point(Generator& generator) {
  again:
    size_t i_region = choose_region(generator);
    Point<F> p = point_generators[i_region](generator);
    for (size_t j = 0; j < i_region; j++) {
      if (region_list[j].shape.contains(p))
        goto again;
    }
    return p;
  }

  template <typename Generator>
  ImageSpaceEventAngle<F> gen_event(Generator& generator) {
    Point<F> p = gen_point(generator);
    F rangle = F(M_PI_4) * uniform_angle(generator);
    return ImageSpaceEventAngle<F>(p.y, p.x, rangle);
  }

  void operator()(int n_emissions) {

    std::vector<std::vector<Event<F>>> event_list_per_thread(
        omp_get_max_threads());

    rng rd;
    std::vector<rng> rng_list;

    for (int i = 0; i < omp_get_max_threads(); ++i) {

      rng_list.push_back(rd);
      rng_list[i].seed(42 + (3453 * i));
      // OR
      // Turn on leapfrogging with an offset that depends on the task id
    }

    for (int i = 0; i < omp_get_max_threads(); ++i)
      event_list_per_thread[i].clear();

#if _OPENMP
#pragma omp for schedule(static)
#endif
    for (int emission = 0; emission < n_emissions; ++emission) {

      auto event = gen_event(rng_list[omp_get_thread_num()]);

      auto res = detector.detect_event(event, rng_list[omp_get_thread_num()]);
      if (res.second) {

        ImageSpaceEventTan<F> revent =
            detector.from_projection_space_tan(res.first);

        if (std::abs(revent.y) >= detector.radius)
          continue;

        Pixel pp = detector.pixel_at(Point<F>(event.z, event.y));

        if (detector.contains_pixel(pp)) {

          Pixel p = detector.pixel_at(Point<F>(revent.z, revent.y));

          output[p.y][p.x]++;
          output_without_errors[pp.y][pp.x]++;

          event_list_per_thread[omp_get_thread_num()].push_back(res.first);
        }
      }
    }

    for (int i = 0; i < omp_get_max_threads(); ++i) {
      events.insert(events.end(),
                    event_list_per_thread[i].begin(),
                    event_list_per_thread[i].end());
    }
  }

  template <class FileWriter>
  void output_bitmap(FileWriter& fw, bool wo_errors = false) {

    fw.template write_header<>(detector.n_z_pixels, detector.n_y_pixels);

    auto& target_output = wo_errors ? output_without_errors : output;
    F output_max = 0;
    for (auto& col : target_output) {
      for (auto& row : col)
        output_max = std::max(output_max, row);
    }

    auto output_gain =
        static_cast<double>(std::numeric_limits<uint8_t>::max()) / output_max;

    uint8_t* row = (uint8_t*)alloca(detector.n_z_pixels);
    for (int y = 0; y < detector.n_y_pixels; ++y) {
      for (auto x = 0; x < detector.n_z_pixels; ++x) {
        row[x] = std::numeric_limits<uint8_t>::max() -
                 output_gain * target_output[y][x];
      }
      fw.write_row(row);
    }
  }

  template <class FileWriter>
  void output_binary(FileWriter& fw, bool wo_errors = false) {

    auto& target_output = wo_errors ? output_without_errors : output;
    for (auto& line : target_output) {
      fw << line;
    }
  }

  template <typename StreamType> Phantom& operator>>(StreamType& out) {
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
}  // Strip
}  // PET2D
