/// Sparse triangle part system matrix binary file format
/// -----------------------------------------------------
/// uint32_t magic = 'PETt'
/// uint32_t n_ pixels_2
/// while (!eof)
///   uint16_t lor_a, lor_b // pair
///   uint32_t pixel_pair_count
///   for(.. count ..)
///     uint16_t t_pixel_x, t_pixel_y
///     uint32_t pixel_hits

#pragma once

#include <map>
#include <random>

#include "detector.h"
#include "circle.h"
#include "bstream.h"
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
  typedef uint32_t file_int;
  typedef uint16_t file_half;
  typedef hit_type *pixels_type;
  typedef pixels_type *matrix_type;
  typedef std::pair<size_t, size_t> lor_type;
  typedef circle<F>   circle_type;
  typedef point<F>    point_type;
  typedef detector<F> detector_type;
  typedef event<F>    event_type;

  /// @param a_n_detectors number of detectors on ring
  /// @param a_n_pixels    number of pixels in each directions
  /// @param a_s_pixel     size of single pixel (pixels are squares)
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
    if(radious    <= 0.)   throw("invalid radious");
    if(w_detector <= 0. ||
       h_detector <= 0.)   throw("invalid detector size");
    if(n_detectors % 4)    throw("number of detectors must be multiple of 4");
    if(n_pixels  % 2)      throw("number of pixels must be multiple of 2");
    if(s_pixel    <= 0.)   throw("invalid pixel size");

    fov_radius=radious/::sqrt(2.0);

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

  /**
   * @param model acceptance model
   *        (returns bool for call operator with given length)
   * @param rx, ry coordinates of the emission point
   * @param output parameter contains the lor of the event
   */
  template <class AcceptanceModel>
  short emit_event(AcceptanceModel &model,
                   F rx, F ry, F angle,
                   lor_type &lor) {

    typename decltype(c_inner)::event_type e(rx, ry, angle);

    // secant for p and phi
    auto i_inner = c_inner.secant_sections(e, n_detectors);
    auto i_outer = c_outer.secant_sections(e, n_detectors);

    auto inner = i_inner.first;
    auto outer = i_outer.first;

    auto hits = 0;

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
      auto step = ((n_detectors+inner-outer) % n_detectors
                   >
                   (n_detectors+outer-inner) % n_detectors
                    ) ? 1 : -1;

#if SKIP_INTERSECTION
      (!(hits++) ? lor.first : lor.second) = i;
#else
      do {
        auto points = (*this)[i].intersections(e);
        // check if we got 2 point intersection
        // then test the model against these points distance

        if ( points.size() == 2 &&
             model( (points[1]-points[0]).length() ) ) {

          hits++;
          (!side ? lor.first : lor.second) = i;
#if COLLECT_INTERSECTIONS
          for(auto &p: points) ipoints.push_back(p);
#endif
          break;
        }
        // step towards outer detector
        prev_i = i, i = (i + step) % n_detectors;

      } while (prev_i != outer); //loop over intersected  detectors
#endif
      if (hits == 0) break;

      // switch side
      inner = i_inner.second;
      outer = i_outer.second;
    }

#if COLLECT_INTERSECTIONS
    if (hits >= 2)
      for(auto &p: ipoints) intersection_points.push_back(p);
#endif

    return hits;
  }

  void add_to_t_matrix(lor_type &lor, size_t i_pixel) {
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

  /// Executes Monte-Carlo system matrix generation for given detector ring
  /// @param gen   random number generator
  /// @param model acceptance model (returns bool for call operator with given length)
  template <class RandomGenerator, class AcceptanceModel>
  void matrix_mc(RandomGenerator &gen, AcceptanceModel model, size_t n_emissions

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

        if ((x*x + y*y) * s_pixel*s_pixel > fov_radius*fov_radius) continue;

        for (auto n = 0; n < n_emissions; ++n) {

          auto rx = ( x + one_dis(gen) ) * s_pixel;
          auto ry = ( y + one_dis(gen) ) * s_pixel;

          // ensure we are within a triangle
          if (rx > ry) continue;

          auto angle= phi_dis(gen);
          lor_type lor;
          auto hits = emit_event(model, rx, ry, angle, lor);

          // do we have hit on both sides?
          if (hits >= 2) {
            auto i_pixel = t_pixel_index(x, y);

            if (o_collect_mc_matrix) {
              add_to_t_matrix(lor,i_pixel);
            }

            if (o_collect_pixel_stats) {
              ++t_hits[i_pixel];
            }

#if COLLECT_INTERSECTIONS
            for(auto &p: ipoints) intersection_points.push_back(p);
#endif
          } //if (hits>=2)
        } // loop over emmisions from pixel
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
    svg << dr.c_outer;
    svg << dr.c_inner;

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

  constexpr static file_int fourcc(int a, int b, int c, int d) {
    return ( (file_int) (((d)<<24) | ((c)<<16) | ((b)<<8) | (a)) );
  };

  // serialization
  static const file_int magic = fourcc('P','E','T','t');

  friend obstream & operator << (obstream &out, detector_ring &dr) {
    out << magic;
    out << static_cast<file_int>(dr.n_pixels_2);

    for (file_half a = 0; a < dr.n_detectors; ++a) {
      for (file_half b = 0; b <= a; ++b) {
        lor_type lor(a, b);
        auto pixels = dr.t_matrix[lor_index(lor)];
        if (pixels) {
          out << a << b;
          // find out count of non-zero pixels
          file_int count = 0;
          for (auto i = 0; i < dr.n_t_matrix_pixels; ++i) {
            if (pixels[i]) count++;
          }
          out << count;

          // write non-zero pixel pairs
          for (file_half y = 0; y < dr.n_pixels_2; ++y) {
            for (file_half x = 0; x <= y; ++x) {
              file_int hits = pixels[t_pixel_index(x, y)];
              if (hits) {
                out << x << y << hits;
              }
            }
          }
        }
      }
    }
    return out;
  }

  friend ibstream & operator >> (ibstream &in, detector_ring &dr) {
    file_int in_magic;
    in >> in_magic;
    if (in_magic != magic) {
      throw("invalid file type format");
    }
    file_int in_n_pixel_2;
    in >> in_n_pixel_2;
    if (in_n_pixel_2 != dr.n_pixels_2) {
      throw("incompatible input matrix dimensions");
    }
    while (!in.eof()) {
      file_half a, b;
      in >> a >> b;
      lor_type lor(a, b);
      auto i_lor = lor_index(lor);
      auto pixels = dr.t_matrix[i_lor];
      if (!pixels) {
        dr.t_matrix[i_lor] = pixels = new hit_type[dr.n_t_matrix_pixels];
      }
      file_int count;
      in >> count;
      // increment hits
      for (auto i = 0; i < count; ++i) {
        file_half x, y;
        file_int hits;
        in >> x >> y >> hits;
        auto i_p = t_pixel_index(x, y);
        pixels[i_p] += hits;
        dr.t_hits[i_p] += hits;
      }
    }
    return in;
  }

private:
  static size_t lor_index(lor_type &lor) {
    if (lor.first < lor.second) {
      std::swap(lor.first, lor.second);
    }
    return lor.first*(lor.first+1)/2 + lor.second;
  }

  static constexpr size_t t_pixel_index(size_t x, size_t y) {
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
  F fov_radius;
};
