/// Sparse system matrix binary file format
/// -----------------------------------------------------
/// uint32_t magic       // 'PETp' (triangular) / 'PETP' (full)
/// uint32_t n_pixels    // half size [PETp] / full size [PETP]
/// uint32_t n_emissions // per pixel
/// uint32_t n_detectors // regardless of magic
/// while (!eof)
///   uint16_t lor_a, lor_b // pair
///   uint32_t pixel_pair_count
///   for(.. count ..)
///     uint16_t pixel_x, pixel_y // half pixels [PETp] / pixels [PETP]
///     uint32_t pixel_hits

#pragma once

#include "detector_ring.h"
#include "bstream.h"

#define fourcc(a, b, c, d) (((d)<<24) | ((c)<<16) | ((b)<<8) | (a))

#if _OPENMP
#include <omp.h>
#endif

/// Provides storage for 1/8 PET system matrix
template <typename F = double, typename HitType = int>
class matrix_mc {
public:
  typedef uint8_t bitmap_pixel_type;
  typedef HitType hit_type;
  typedef uint32_t file_int;
  typedef uint16_t file_half;
  typedef hit_type *pixels_type;
  typedef pixels_type *matrix_type;
  typedef typename detector_ring<F>::lor_type lor_type;

  /// @param a_dr        detector ring (model)
  /// @param a_n_pixels  number of pixels in each directions
  /// @param a_s_pixel   size of single pixel (pixels are squares)
  matrix_mc(detector_ring<F> &a_dr, size_t a_n_pixels, F a_s_pixel)
    : dr(a_dr)
    , n_pixels(a_n_pixels)
    , n_pixels_2(a_n_pixels/2)
    , n_t_matrix_pixels( a_n_pixels/2 * (a_n_pixels/2+1) / 2 )
    , s_pixel(a_s_pixel)
    , n_emissions(0)
    , n_2_detectors(2 * dr.detectors())
    , n_1_detectors_2(dr.detectors() + dr.detectors() / 2)
    , n_1_detectors_4(dr.detectors() + dr.detectors() / 4)
    , n_lors( dr.lors() )
    , output_triangular(true) {
    if(n_pixels % 2 ) throw("number of pixels must be multiple of 2");
    if(s_pixel <= 0.) throw("invalid pixel size");

    // reserve for all lors
    t_matrix = new pixels_type[n_lors]();

    // reserve for pixel stats
    t_hits = new hit_type[n_t_matrix_pixels]();
  }

  ~matrix_mc() {
    for (auto i = 0; i < n_lors; ++i) {
      if (t_matrix[i]) delete [] t_matrix[i];
    }
    delete [] t_matrix;
    delete [] t_hits;
  }

  F     pixel_size()    const {return s_pixel;}
  int   get_n_pixels()  const {return n_pixels;}

  void add_to_t_matrix(lor_type &lor, size_t i_pixel) {
    auto i_lor  = t_lor_index(lor);
    auto pixels = t_matrix[i_lor];

    // prealocate pixels for specific lor
    // all are nulls upon startup
    if (!pixels) {
#if _OPENMP
      // entering critical section, and check again
      // because it may have changed in meantime
#pragma omp critical
      if ( !(pixels = t_matrix[i_lor]) ) {
        pixels = t_matrix[i_lor] = new hit_type[n_t_matrix_pixels]();
      }
#else
      pixels = t_matrix[i_lor] = new hit_type[n_t_matrix_pixels]();
#endif
    }
    ++pixels[i_pixel];
  }

  hit_type operator () (lor_type lor, size_t x, size_t y) const {
    bool diag; int symmetry;
    auto i_pixel = pixel_index(x, y, diag, symmetry);
    auto pixels = t_matrix[lor_index(lor, symmetry)];
    return pixels ? pixels[i_pixel] * (diag ? 2 : 1) : 0;
  }

  void add_hit(int i_pixel) {
    ++t_hits[i_pixel];
  }

  F non_zero_lors() {
    size_t non_zero_lors = 0;
    for (auto i = 0; i < n_lors; ++i) {
      if (t_matrix[i]) ++non_zero_lors;
    }
    return non_zero_lors;
  }

  hit_type hits(size_t x, size_t y) const {
    bool diag; int symmetry;
    auto i_pixel = pixel_index(x, y, diag, symmetry);
    return t_hits[i_pixel] * (diag ? 2 : 1);
  }

  void increase_n_emissions(int  count) {n_emissions+=count;}

  static constexpr size_t t_pixel_index(size_t x, size_t y) {
    return y*(y+1)/2 + x;
  }

  /**Computes pixel index and determines symmetry number based on pixel
     position
     @param x        pixel x coordinate (0..n_pixels)
     @param x        pixel y coordinate (0..n_pixels)
     @param diag     outputs true if abs(x)==abs(y)
     @param symmetry outputs symmetry number (0..7)
  */
  size_t pixel_index(ssize_t x, ssize_t y, bool &diag, int &symmetry) const {
    // shift so 0,0 is now center
    x -= n_pixels_2; y -= n_pixels_2;
    // mirror
    symmetry = 0;
    if (x < 0) {
      x = -x-1;
      symmetry |= 2;
    };
    if (y < 0) {
      y = -y-1;
      symmetry |= 1;
    }
    // triangulate
    if (x > y) {
      std::swap(x, y);
      symmetry |= 4;
    };
    diag = (x == y);
    return t_pixel_index(x, y);
  }


  template<class FileWriter>
  void output_bitmap(FileWriter &fw) {
    fw.template write_header<bitmap_pixel_type>(n_pixels, n_pixels);
    hit_type pixel_max = 0;
    for (auto y = 0; y < n_pixels; ++y) {
      for (auto x = 0; x < n_pixels; ++x) {
        pixel_max = std::max(pixel_max, hits(x, y));
      }
    }
    auto gain = static_cast<double>(std::numeric_limits<bitmap_pixel_type>::max()) / pixel_max;
    for (auto y = 0; y < n_pixels; ++y) {
      bitmap_pixel_type row[n_pixels];
      for (auto x = 0; x < n_pixels; ++x) {
        row[x] = std::numeric_limits<bitmap_pixel_type>::max() - gain * hits(x, y);
      }
      fw.write_row(row);
    }
  }

  template<class FileWriter>
  void output_bitmap(FileWriter &fw, lor_type &lor) {
    fw.template write_header<bitmap_pixel_type>(n_pixels, n_pixels);
    hit_type pixel_max = 0;
    for (auto y = 0; y < n_pixels; ++y) {
      for (auto x = 0; x < n_pixels; ++x) {
        pixel_max = std::max(pixel_max, (*this)(lor, x, y));
      }
    }
    auto gain = static_cast<double>(std::numeric_limits<bitmap_pixel_type>::max()) / pixel_max;
    for (int y = n_pixels-1; y >= 0; --y) {
      bitmap_pixel_type row[n_pixels];
      for (auto x = 0; x < n_pixels; ++x) {
        row[x] = std::numeric_limits<bitmap_pixel_type>::max() - gain * (*this)(lor, x, y);
      }
      fw.write_row(row);
    }
  }

  // binary serialization                                  // n_pixels  n_detectors  triagular
  static const file_int magic_1 = fourcc('P','E','T','t'); //                           X
  static const file_int magic_2 = fourcc('P','E','T','s'); //     X                     X
  static const file_int magic_t = fourcc('P','E','T','p'); //     X          X          X
  static const file_int magic_f = fourcc('P','E','T','P'); //     X          X

  friend obstream & operator << (obstream &out, matrix_mc &mmc) {
    if (mmc.output_triangular) {
      out << magic_t;
      out << static_cast<file_int>(mmc.n_pixels_2);
    } else {
      out << magic_f;
      out << static_cast<file_int>(mmc.n_pixels);
    }
    out << static_cast<file_int>(mmc.n_emissions);
    out << static_cast<file_int>(mmc.dr.detectors());

    for (file_half a = 0; a < mmc.dr.detectors(); ++a) {
      for (file_half b = 0; b <= a; ++b) {
        lor_type lor(a, b);
        if (mmc.output_triangular) {
          auto pixels = mmc.t_matrix[t_lor_index(lor)];
          if (pixels) {
            out << a << b;
            // find out count of non-zero pixels
            file_int count = 0;
            for (auto i = 0; i < mmc.n_t_matrix_pixels; ++i) {
              if (pixels[i]) count++;
            }
            out << count;

            // write non-zero pixel pairs
            for (file_half y = 0; y < mmc.n_pixels_2; ++y) {
              for (file_half x = 0; x <= y; ++x) {
                file_int hits = pixels[t_pixel_index(x, y)];
                if (hits) {
                  out << x << y << hits;
                }
              }
            }
          }
        } else { // output full (may be slow)
          // find out count of non-zero pixels
          file_int count = 0;
          for (file_half x = 0; x < mmc.n_pixels; ++x) {
            for (file_half y = 0; y < mmc.n_pixels; ++y) {
              if (mmc(lor, x, y)) count++;
            }
          }
          if (count) {
            out << a << b << count;
            // write out non-zero hits
            for (file_half x = 0; x < mmc.n_pixels; ++x) {
              for (file_half y = 0; y < mmc.n_pixels; ++y) {
                file_int hits = mmc(lor, x, y);
                if (hits) {
                  out << x << y << hits;
                }
              }
            }
          }
        }
      }
    }
    return out;
  }

  friend ibstream & operator >> (ibstream &in, matrix_mc &mmc) {
    // read header
    file_int in_is_triangular, in_n_pixels, in_n_emissions, in_n_detectors;
    mmc.read_header(in, in_is_triangular, in_n_pixels, in_n_emissions, in_n_detectors);

    if (mmc.n_emissions && !in_is_triangular) {
      throw("full matrix cannot be loaded to non-empty matrix");
    }
 
    // validate incoming parameters
    if (in_n_pixels && in_n_pixels != (in_is_triangular ? mmc.n_pixels_2 : mmc.n_pixels)) {
      std::ostringstream msg;
      msg << "incompatible input matrix dimensions "
          << in_n_pixels
          << " != "
          << mmc.n_pixels_2;
      throw(msg.str());
    }
    if (in_n_detectors && in_n_detectors != mmc.dr.detectors()) {
      throw("incompatible input number of detectors");
    }

    mmc.n_emissions += in_n_emissions;

    // load hits
    for (;;) {
      file_half a, b;
      in >> a >> b;
      lor_type lor(a, b);

      if (in.eof()) break;

      file_int count;
      in >> count;

      if (in_is_triangular) {
        auto i_lor = t_lor_index(lor);
        if (i_lor >= mmc.n_lors) {
          std::ostringstream msg;
          msg << "invalid LOR address ("
              << a << "," << b << ")";
          throw(msg.str());
        }

        auto pixels = mmc.t_matrix[i_lor];
        if (!pixels) {
          mmc.t_matrix[i_lor] = pixels = new hit_type[mmc.n_t_matrix_pixels]();
        }
        // increment hits
        for (auto i = 0; i < count; ++i) {
          file_half x, y;
          file_int hits;
          in >> x >> y >> hits;
          auto i_pixel = t_pixel_index(x, y);
          pixels[i_pixel]    += hits;
          mmc.t_hits[i_pixel] += hits;
        }
      } else { // full matrix
        // increment hits
        for (auto i = 0; i < count; ++i) {
          file_half x, y;
          file_int hits;
          in >> x >> y >> hits;

          bool diag; int symmetry;
          auto i_pixel = mmc.pixel_index(x, y, diag, symmetry);
          if (i_pixel >= mmc.n_t_matrix_pixels) throw("invalid pixel address");

          auto i_lor   = mmc.lor_index(lor, symmetry);
          if (i_lor >= mmc.n_lors) {
            std::ostringstream msg;
            msg << "invalid LOR address ("
                << a << "," << b;
            throw(msg.str());
          }

          auto pixels  = mmc.t_matrix[i_lor];
          if (!pixels) {
            mmc.t_matrix[i_lor] = pixels = new hit_type[mmc.n_t_matrix_pixels]();
          }
          pixels[i_pixel]    += hits;
          mmc.t_hits[i_pixel] += hits;
        }
      }
    }

    // we need to divide all values by 8 when loading full matrix
    if (!in_is_triangular) {
      for (auto i_lor = 0; i_lor < mmc.n_lors; ++i_lor) {
        auto pixels  = mmc.t_matrix[i_lor];
        if (!pixels) continue;
        for (auto i_pixel = 0; i_pixel < mmc.n_t_matrix_pixels; ++i_pixel) {
          pixels[i_pixel] /= 8;
        }
      }
    }

    return in;
  }

  // text output (for validation)
  friend std::ostream & operator << (std::ostream &out, matrix_mc &mmc) {
    out << "n_pixels="    << mmc.n_pixels    << std::endl;
    out << "n_emissions=" << mmc.n_emissions << std::endl;
    out << "n_detectors=" << mmc.dr.detectors() << std::endl;

    for (file_half a = 0; a < mmc.dr.detectors(); ++a) {
      for (file_half b = 0; b <= a; ++b) {
        lor_type lor(a, b);
        auto pixels = mmc.t_matrix[t_lor_index(lor)];
        if (pixels) {
          out << "  lor=(" << a << "," << b << ")" << std::endl;
          // find out count of non-zero pixels
          file_int count = 0;
          for (auto i = 0; i < mmc.n_t_matrix_pixels; ++i) {
            if (pixels[i]) count++;
          }
          out << "  count=" << count << std::endl;

          // write non-zero pixel pairs
          for (file_half y = 0; y < mmc.n_pixels_2; ++y) {
            for (file_half x = 0; x <= y; ++x) {
              file_int hits = pixels[t_pixel_index(x, y)];
              if (hits) {
                out << "    (" << x << "," << y << ")=" << hits << std::endl;
              }
            }
          }
        }
      }
    }
    return out;
  }

  static void read_header(ibstream &in
                          , file_int &in_is_triangular
                          , file_int &in_n_pixels
                          , file_int &in_n_emissions
                          , file_int &in_n_detectors
                          ) {
    file_int in_magic;
    in >> in_magic;
    if (in_magic != magic_t && in_magic != magic_f && in_magic != magic_1 && in_magic != magic_2) {
      throw("invalid file type format");
    }
    in_is_triangular = (in_magic != magic_f);

    // load matrix size
    in >> in_n_pixels;

    // load number of emissions
    if (in_magic == magic_t || in_magic == magic_f || in_magic == magic_2) {
      in >> in_n_emissions;
    }

    // load number of detectors
    in_n_detectors = 0;
    if (in_magic == magic_t || in_magic == magic_f) {
      in >> in_n_detectors;
    }
  }

  bool output_triangular;

private:
  static size_t t_lor_index(lor_type &lor) {
    if (lor.first < lor.second) {
      std::swap(lor.first, lor.second);
    }
    return lor.first*(lor.first+1)/2 + lor.second;
  }

  /// Computes LOR index based on given symmetry (1 out 8)
  /// @param lor      detector number pair
  /// @param symmetry number (0..7)
  size_t lor_index(lor_type lor, int symmetry) const {
    if (symmetry & 1) {
      lor.first  = ( n_2_detectors - lor.first  ) % dr.detectors();
      lor.second = ( n_2_detectors - lor.second ) % dr.detectors();
    }
    if (symmetry & 2) {
      lor.first  = ( n_1_detectors_2 - lor.first  ) % dr.detectors();
      lor.second = ( n_1_detectors_2 - lor.second ) % dr.detectors();
    }
    if (symmetry & 4) {
      lor.first  = ( n_1_detectors_4 - lor.first  ) % dr.detectors();
      lor.second = ( n_1_detectors_4 - lor.second ) % dr.detectors();
    }
    return t_lor_index(lor);
  }


  detector_ring<F> &dr;
  matrix_type t_matrix;
  pixels_type t_hits;
  size_t n_t_matrix_pixels;
  hit_type n_emissions;
  size_t n_pixels;
  size_t n_pixels_2;
  size_t n_2_detectors;
  size_t n_1_detectors_2;
  size_t n_1_detectors_4;
  size_t n_lors;
  F s_pixel;
};
