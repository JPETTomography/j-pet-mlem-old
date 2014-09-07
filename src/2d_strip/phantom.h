#pragma once

#include <cmath>
#include <vector>
#include <ctime>
#include <random>
#include <algorithm>

#include "event.h"
#include "geometry/ellipse.h"

#if _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#endif

template <typename F> int sgn(F val) { return (0 < val) - (val < 0); }

template <typename F> class PhantomRegion {
 public:
  PhantomRegion(const Ellipse<F>& ellipse, F intensity)
      : ellipse_(ellipse), intensity_(intensity) {
    weight_ = intensity_ * ellipse_.measure();
  }

  bool in(F x, F y) const { return ellipse_.in(x, y); }
  F weight() const { return weight_; }
  Ellipse<F> shape() const { return ellipse_; }

 private:
  Ellipse<F> ellipse_;
  F intensity_;
  F weight_;
};

template <typename D, typename FType = double> class Phantom {
  typedef FType F;
  typedef ::Pixel<> Pixel;
  typedef std::minstd_rand0 rng;

 private:
  D detector_;

  std::vector<PhantomRegion<F>> region_list;
  std::vector<F> CDF_;
  std::vector<EllipsePointsGenerator<F>> point_generators_;

  std::vector<Event<F>> events;
  std::vector<std::vector<F>> output;
  std::vector<std::vector<F>> output_without_errors;

  std::uniform_real_distribution<F> uni_;
  std::uniform_real_distribution<F> uniform_angle;

 public:
  Phantom(const D& detector, const std::vector<PhantomRegion<F>>& el)
      : detector_(detector), CDF_(el.size(), 0), uniform_angle(-1, 1) {
    region_list = el;
    CDF_[0] = region_list[0].weight();

    for (size_t i = 1; i < el.size(); i++) {
      CDF_[i] = region_list[i].weight() + CDF_[i - 1];
    }
    F norm = CDF_[el.size() - 1];
    for (size_t i = 0; i < el.size(); i++) {
      CDF_[i] /= norm;
    }

    for (size_t i = 0; i < el.size(); ++i)
      point_generators_.push_back(EllipsePointsGenerator<F>(el[i].shape()));

    output.assign(detector_.n_y_pixels, std::vector<F>(detector_.n_z_pixels, 0));
    output_without_errors.assign(detector_.n_y_pixels,
                                 std::vector<F>(detector_.n_z_pixels, 0));
  }

  template <typename G> size_t choose_region(G& gen) {
    F r = uni_(gen);
    size_t i = 0;

    while (r > CDF_[i])
      ++i;

    return i;
  }

  template <typename G> Point<F> gen_point(G& gen) {
  again:
    size_t i_region = choose_region(gen);
    Point<F> p = point_generators_[i_region].point(gen);
    for (int j = 0; j < i_region; j++) {
      if (region_list[j].shape().in(p.x, p.y))
        goto again;
    }
    return p;
  }

  template <typename G> ImageSpaceEventAngle<F> gen_event(G& gen) {
    Point<F> p = gen_point(gen);
    F rangle = F(M_PI_4) * uniform_angle(gen);
    return ImageSpaceEventAngle<F>(p.y, p.x, rangle);
  }

  void operator()(size_t n_emissions) {

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

      auto res = detector_.detect_event(event, rng_list[omp_get_thread_num()]);
      if (res.second) {

        ImageSpaceEventTan<F> revent =
            detector_.from_projection_space_tan(res.first);

        Pixel pp = detector_.pixel_location(Point<F>(event.z, event.y));
        Pixel p = detector_.pixel_location(Point<F>(revent.z, revent.y));

        output[p.y][p.x]++;
        output_without_errors[pp.y][pp.x]++;

        event_list_per_thread[omp_get_thread_num()].push_back(res.first);
      }
    }

    for (int i = 0; i < omp_get_max_threads(); ++i) {
      events.insert(events.end(),
                    event_list_per_thread[i].begin(),
                    event_list_per_thread[i].end());
    }

    std::cout << "Detected  " << events.size() << " events" << std::endl;
  }

  template <class FileWriter>
  void output_bitmap(FileWriter& fw, bool wo_errors = false) {

    fw.template write_header<>(detector_.n_z_pixels, detector_.n_y_pixels);

    auto& target_output = wo_errors ? output_without_errors : output;
    F output_max = 0;
    for (auto& col : target_output) {
      for (auto& row : col)
        output_max = std::max(output_max, row);
    }

    auto output_gain =
        static_cast<double>(std::numeric_limits<uint8_t>::max()) / output_max;

    uint8_t* row = (uint8_t*)alloca(detector_.n_z_pixels);
    for (int y = 0; y < detector_.n_y_pixels; ++y) {
      for (auto x = 0; x < detector_.n_z_pixels; ++x) {
        row[x] = std::numeric_limits<uint8_t>::max() -
                 output_gain * target_output[y][x];
      }
      fw.write_row(row);
    }
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
