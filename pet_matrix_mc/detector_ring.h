#pragma once

#include <map>
#include <random>

#include "detector.h"
#include "circle.h"
#include "svg_ostream.h"

#if _OPENMP
#include <omp.h>
#endif

/// Provides model for 2D ring of detectors
template <typename F = double, typename HitType = int>
class detector_ring : public std::vector<detector<F>> {
public:
  typedef uint8_t bitmap_pixel_type;
  typedef HitType hit_type;
  typedef hit_type *pixels_type;
  typedef pixels_type *matrix_type;
  typedef std::pair<size_t, size_t> lor_type;
  typedef circle<F> circle_type;
  typedef point<F> point_type;
  typedef detector<F> detector_type;

  /// @param a_n_detectors number of detectors on ring
  /// @param a_n_pixels    number of pixels in directions
  /// @param a_s_pixel     size of single pixel
  /// @param radious       radious of ring
  /// @param w_detector    width of single detector (along ring)
  /// @param h_detector    height/depth of single detector (perpendicular to ring)
  detector_ring(size_t a_n_detectors, size_t a_n_pixels, F a_s_pixel, F radious, F w_detector, F h_detector)
  : c_inner(radious)
  , c_outer(radious+h_detector)
  , n_detectors(a_n_detectors)
  , n_pixels(a_n_pixels)
  , n_pixels_2(a_n_pixels/2)
  , n_t_matrix_pixels( a_n_pixels/2 * (a_n_pixels/2+1) / 2 )
  , n_lors( a_n_detectors * (a_n_detectors+1) / 2 )
  , s_pixel(a_s_pixel)
  {
    // reserve for all lors
    t_matrix = new pixels_type[n_lors]();

    // reserve for pixel stats
    t_hits = new hit_type[n_t_matrix_pixels];

    detector_type detector_base(h_detector, w_detector);

    // move detector to the right edge of inner ring
    // along zero angle polar coordinate
    detector_base += point<>(radious + h_detector/2, 0.);

    // produce detector ring rotating base detector n times
    std::vector<detector<>> detector_ring;
    for (auto n = 0; n < n_detectors; ++n) {
      this->push_back( detector_base.rotated(2. * M_PI * n / n_detectors) );
    }
  }

  ~detector_ring() {
    for (auto i = 0; i < n_lors; ++i) {
      if (t_matrix[i]) delete [] t_matrix[i];
    }
    delete [] t_matrix;
    delete [] t_hits;
  }

  /// Executes Monte-Carlo system matrix generation for given detector ring
  /// @param gen   random number generator
  /// @param model acceptance model (returns bool for call operator with given length)
  template <class RandomGenerator, class AcceptanceModel>
  void matrix_mc(RandomGenerator gen, AcceptanceModel model, size_t n_emissions

  , bool o_collect_mc_matrix = true
  , bool o_collect_pixel_stats = true) {

    std::uniform_real_distribution<> one_dis(0., 1.);
    std::uniform_real_distribution<> phi_dis(0., M_PI);
    // iterating only triangular matrix,
    // being upper right part or whole system matrix
#if _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    // descending, since biggest chunks start first, but may end last
    for (ssize_t y = n_pixels_2 - 1; y >= 0; --y) {
      for (auto x = 0; x <= y; ++x) {
        // if (x > y) continue;
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
          auto inner = i_inner.first;
          auto outer = i_inner.first;

          auto hits = 0;
          lor_type lor;
#if COLLECT_INTERSECTIONS
          std::vector<point_type> ipoints;
#endif
          // process possible hit detectors on both sides
          for (auto side = 0; side < 2; ++side) {
            // starting from inner index
            // iterating to outer index
            auto i = inner;
            auto prev_i = i;

            // tells in which direction we got shorter modulo distance
            auto step = (n_detectors+inner-outer) % n_detectors
                      > (n_detectors+outer-inner) % n_detectors ? 1 : -1;

#if SKIP_INTERSECTION
            (!(hits++) ? lor.first : lor.second) = i;
#else
            do {
              auto points = (*this)[i].intersections(e);
              // check if we got 2 point intersection
              // then test the model against these points distance
              if ( points.size() == 2 &&
                   model( (points[1]-points[0]).length() ) ) {
                (!(hits++) ? lor.first : lor.second) = i;
#if COLLECT_INTERSECTIONS
                for(auto &p: points) ipoints.push_back(p);
#endif
                break;
              }
              // step towards outer detector
              prev_i = i, i = (i + step) % n_detectors;
            } while (prev_i != outer);
#endif
            // switch side
            inner = i_inner.second;
            outer = i_outer.second;
          }

          // do we have hit on both sides?
          if (hits >= 2) {
            auto i_pixel = t_pixel_index(x, y);

            if (o_collect_mc_matrix) {
              auto i_lor  = lor_index(lor);
              auto pixels = t_matrix[i_lor];

              // prealocate pixels for specific lor
              // all are nulls upon startup
              if (!pixels) {
#if _OPENMP
                // entering critical section, and check again
                // because it may have changed in meantime
                #pragma omp critical
                if ( !(pixels = t_matrix[i_lor]) ) {
                  pixels = t_matrix[i_lor] = new hit_type[n_t_matrix_pixels];
                }
#else
                t_matrix[i_lor] = pixels = new hit_type[n_t_matrix_pixels];
#endif
              }
              ++pixels[i_pixel];
            }

            if (o_collect_pixel_stats) {
              ++t_hits[i_pixel];
            }

#if COLLECT_INTERSECTIONS
            for(auto &p: ipoints) intersection_points.push_back(p);
#endif
          }
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

  hit_type matrix(lor_type lor, size_t x, size_t y) const {
    auto pixels = t_matrix[lor_index(lor)];
    bool diag;
    return pixels ? pixels[pixel_index(x, y, diag)] * (diag ? 2 : 1) : 0;
  }

  hit_type matrix(size_t lor, size_t x, size_t y) const {
    auto pixels = t_matrix[lor];
    bool diag;
    return pixels ? pixels[pixel_index(x, y, diag)] * (diag ? 2 : 1) : 0;
  }

  hit_type hits(size_t x, size_t y) const {
    bool diag;
    return t_hits[pixel_index(x, y, diag)] * (diag ? 2 : 1);
  }

  template<class FileWriter>
  void output_bitmap(FileWriter &fw, hit_type pixel_max, ssize_t lor) {
    fw.template write_header<bitmap_pixel_type>(n_pixels, n_pixels);
    auto gain = static_cast<double>(std::numeric_limits<bitmap_pixel_type>::max()) / pixel_max;
    for (auto y = 0; y < n_pixels; ++y) {
      bitmap_pixel_type row[n_pixels];
      for (auto x = 0; x < n_pixels; ++x) {
        auto v = (lor > 0 ? matrix(lor, x, y) : hits(x, y));
        row[x] = std::numeric_limits<bitmap_pixel_type>::max() - gain * v;
      }
      fw.write_row(row);
    }
  }

  friend svg_ostream<F> & operator << (svg_ostream<F> &svg, detector_ring &dr) {
    svg << dr.c_inner;
    svg << dr.c_outer;

    for (auto detector: dr) {
      svg << detector;
    }
#if COLLECT_INTERSECTIONS
    for (auto &p: dr.intersection_points) {
      svg << "<circle cx=\"" << p.x << "\" cy=\"" << p.y << "\" r=\"0.002\"/>" << std::endl;
    }
#endif
    return svg;
  }

private:
  size_t lor_index(lor_type &lor) const {
    if (lor.first < lor.second) {
      std::swap(lor.first, lor.second);
    }
    return lor.first*(lor.first+1)/2 + lor.second;
  }

  size_t t_pixel_index(size_t x, size_t y) const {
    return y*(y+1)/2 + x;
  }

  size_t pixel_index(ssize_t x, ssize_t y, bool &diag) const {
    // shift so 0,0 is now center
    x -= n_pixels_2; y -= n_pixels_2;
    // mirror
    if (x < 0) x = -x-1;
    if (y < 0) y = -y-1;
    // triangulate
    if (x > y) std::swap(x, y);
    diag = (x == y);
    return t_pixel_index(x, y);
  }

private:
#if COLLECT_INTERSECTIONS
  std::vector<point_type> intersection_points;
#endif
  matrix_type t_matrix;
  pixels_type t_hits;
  size_t n_t_matrix_pixels;
  circle_type c_inner;
  circle_type c_outer;
  size_t n_pixels;
  size_t n_pixels_2;
  size_t n_detectors;
  size_t n_lors;
  F s_pixel;
};
