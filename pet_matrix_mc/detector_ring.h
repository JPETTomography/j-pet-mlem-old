#pragma once

#include <map>
#include <random>

#include "detector.h"
#include "circle.h"

template <typename F = double, typename HitType = int>
class detector_ring : public std::vector<detector<F>> {
public:
  typedef HitType hit_type;
  typedef hit_type *pixels_type;
  typedef pixels_type *matrix_type;
  typedef std::pair<size_t, size_t> lor_type;
  typedef circle<F> circle_type;
  typedef detector<F> detector_type;

  detector_ring(size_t a_n_detectors, size_t a_n_pixels, F a_s_pixel, F radious, F w_detector, F h_detector)
  : c_inner(radious)
  , c_outer(radious+h_detector)
  , n_detectors(a_n_detectors)
  , n_pixels(a_n_pixels)
  , n_t_matrix_pixels( a_n_pixels/2 * (a_n_pixels/2+1) / 2 )
  , n_lors( a_n_detectors * (a_n_detectors+1) / 2 )
  , s_pixel(a_s_pixel)
  {
    // reserve for all lors
    t_matrix = new pixels_type[n_lors]();

    detector_type detector_base(h_detector, w_detector);

    // move detector to the edge of inner ring
    detector_base += point<>(radious + h_detector/2, 0.);

    // produce detector ring
    std::vector<detector<>> detector_ring;
    for (auto n = 0; n < n_detectors; ++n) {
      this->push_back( detector_base.rotated(2. * M_PI / n_detectors) );
    }
  }

  ~detector_ring() {
    for (auto i = 0; i < n_lors; ++i) {
      if (t_matrix[i]) delete [] t_matrix[i];
    }
    delete [] t_matrix;
  }

  template <class RandomGenerator, class AcceptanceModel>
  void matrix_mc(RandomGenerator gen, AcceptanceModel model, size_t n_emissions) {
    std::uniform_real_distribution<> one_dis(0, 1);
    std::uniform_real_distribution<> phi_dis(0, M_PI);

    // iterating only triangular matrix,
    // being upper right part or whole system matrix
    for (auto y = 0; y < n_pixels/2; ++y)
      for (auto x = 0; x <= y; ++x)
        for (auto n = 0; n < n_emissions; ++n) {
          auto rx = ( x + one_dis(gen) ) * s_pixel;
          auto ry = ( y + one_dis(gen) ) * s_pixel;
          // ensure we are within a triangle
          if (rx > ry) continue;
          // random point within a pixel
          typename decltype(c_inner)::event_type e(rx, ry, phi_dis(gen));
          // secant for p and phi
          auto i_inner = c_inner.secant_sections(e, n_detectors);
          auto i_outer = c_outer.secant_sections(e, n_detectors);
#if DEBUG
          std::cout << '(' << e.x << ',' << e.y << '@' << e.phi << ')' << std::endl;
#endif
          auto inner = i_inner.first;
          auto outer = i_inner.first;

          auto hits = 0;
          lor_type lor;

          // process both sides
          for (auto side = 0; side < 2; ++side) {
            // starting from inner index
            // iterating to outer index
            auto i = inner;
            auto prev_i = i;
            auto step = inner <= outer ? 1 : -1;
#if DEBUG
            std::cout << '[' << side << "] " << inner  << '-' << outer << ':' << step;
#endif
            do {
              auto points = (*this)[i].intersections(e);
              // bail out if intersects
              if ( points.size() == 2 &&
                   model( (points[1]-points[0]).length() ) ) {
                // FIXME: we shall count probability here now!
#if DEBUG
                std::cout << ' ' << i;
#endif
                (!(hits++) ? lor.first : lor.second) = i;
                break;
              }
              // step
              prev_i = i, i = (i + step) % n_detectors;
#if DEBUG
              std::cout << " s<" << i;
#endif
            } while (prev_i != outer);

            // switch side
            inner = i_inner.second;
            outer = i_outer.second;
#if DEBUG
            std::cout << std::endl;
#endif
            if (hits >= 2) {
              lor = std::make_pair(
                std::max(lor.first, lor.second),
                std::min(lor.first, lor.second)
              );
              auto i_lor   = lor.first*(lor.first+1)/2 + lor.second;
              auto i_pixel = y*(y+1)/2 + x;
              auto pixels  = t_matrix[i_lor];

              // prealocate pixels for specific lor
              if (!pixels) {
                t_matrix[i_lor] = pixels = new hit_type[n_t_matrix_pixels];
              }

              ++pixels[i_pixel];
            }
          }
        }
  }

  size_t lors() const { return n_lors; }

  F non_zero_lors() {
    size_t non_zero_lors = 0;
    for (auto i = 0; i < n_lors; ++i) {
      if (t_matrix[i]) ++non_zero_lors;
    }
    return non_zero_lors;
  }

private:
  matrix_type t_matrix;
  size_t n_t_matrix_pixels;
  circle_type c_inner;
  circle_type c_outer;
  size_t n_pixels;
  size_t n_detectors;
  size_t n_lors;
  F s_pixel;
};
