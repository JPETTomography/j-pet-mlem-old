/// Sparse system matrix binary file format
/// -----------------------------------------------------
/// uint32_t magic       // 'PETp' (triangular) / 'PETP' (full)
/// uint32_t n_pixels_    // half size [PETp] / full size [PETP]
/// uint32_t n_emissions // per pixel
/// uint32_t n_detectors // regardless of magic
/// while (!eof)
///   uint16_t lor_a, lor_b // pair
///   uint32_t pixel_pair_count
///   for(.. count ..)
///     uint16_t pixel_x, pixel_y // half pixels [PETp] / pixels [PETP]
///     uint32_t pixel_hits

#pragma once

#if _OPENMP
#include <omp.h>
#endif

#include "detector_ring.h"
#include "bstream.h"
#include "triangular_pixel_map.h"

#define fourcc(a, b, c, d) (((d) << 24) | ((c) << 16) | ((b) << 8) | (a))

/// Provides storage for 1/8 PET system matrix
template <typename F = double, typename HitType = int>
class MatrixLORMajor : public TriangularPixelMap<F, HitType> {
 public:
  typedef HitType Hit;
  typedef uint8_t BitmapPixel;
  typedef uint32_t FileInt;
  typedef uint16_t FileHalf;
  typedef Hit* Pixels;
  typedef Pixels* Matrix;
  typedef typename DetectorRing<F>::LOR LOR;
  typedef TriangularPixelMap<F, Hit> Super;

  /// @param dr       detector ring (model)
  /// @param n_pixels number of pixels in each directions
  /// @param s_pixel  size of single pixel (pixels are squares)
  MatrixLORMajor(DetectorRing<F>& dr, size_t n_pixels, F s_pixel)
      : TriangularPixelMap<F, Hit>(n_pixels),
        output_triangular(true),
        dr_(dr),
        n_emissions_(0),
        n_2_detectors_(2 * dr.detectors()),
        n_1_detectors_2_(dr.detectors() + dr.detectors() / 2),
        n_1_detectors_4_(dr.detectors() + dr.detectors() / 4),
        n_lors_(dr.lors()),
        s_pixel_(s_pixel) {
    if (n_pixels % 2)
      throw("number of pixels must be multiple of 2");
    if (s_pixel <= 0.)
      throw("invalid pixel size");

    // reserve for all lors
    t_matrix_ = new Pixels[n_lors_]();
  }

  ~MatrixLORMajor() {
    for (auto i = 0; i < n_lors_; ++i) {
      if (t_matrix_[i])
        delete[] t_matrix_[i];
    }
    delete[] t_matrix_;
  }

  F pixel_size() const { return s_pixel_; }
  int get_n_pixels() const { return Super::n_pixels_in_row(); }

  void add_to_t_matrix(LOR& lor, size_t i_pixel) {
    auto i_lor = t_lor_index(lor);
    auto pixels = t_matrix_[i_lor];

    // prealocate pixels for specific lor
    // all are nulls upon startup
    if (!pixels) {
#if _OPENMP
// entering critical section, and check again
// because it may have changed in meantime
#pragma omp critical
      if (!(pixels = t_matrix_[i_lor])) {
        pixels = t_matrix_[i_lor] =
                 new Hit[Super::total_n_pixels_in_triangle()]();
      }
#else
      pixels = t_matrix_[i_lor] =
               new Hit[Super::total_n_pixels_in_triangle()]();
#endif
    }
    ++pixels[i_pixel];
  }

  Hit operator()(LOR lor, size_t x, size_t y) const {
    bool diag;
    int symmetry;
    auto i_pixel = this->pixel_index(x, y, diag, symmetry);
    auto pixels = t_matrix_[lor_index(lor, symmetry)];
    return pixels ? pixels[i_pixel] * (diag ? 2 : 1) : 0;
  }

  F non_zero_lors() {
    size_t non_zero_lors = 0;
    for (auto i = 0; i < n_lors_; ++i) {
      if (t_matrix_[i])
        ++non_zero_lors;
    }
    return non_zero_lors;
  }

  void increase_n_emissions(int n_emissions) { n_emissions_ += n_emissions; }

  // binary serialization                 // n_pixels_  n_detectors  triagular
  static const FileInt MAGIC_VERSION_1 =
      fourcc('P', 'E', 'T', 't');  //                           X
  static const FileInt MAGIC_VERSION_2 =
      fourcc('P', 'E', 'T', 's');  //     X                     X
  static const FileInt MAGIC_VERSION_TRIANGULAR =
      fourcc('P', 'E', 'T', 'p');  //     X          X          X
  static const FileInt MAGIC_VERSION_FULL =
      fourcc('P', 'E', 'T', 'P');  //     X          X

  friend obstream& operator<<(obstream& out, MatrixLORMajor& mmc) {
    if (mmc.output_triangular) {
      out << MAGIC_VERSION_TRIANGULAR;
      out << static_cast<FileInt>(mmc.n_pixels_in_row_half());
    } else {
      out << MAGIC_VERSION_FULL;
      out << static_cast<FileInt>(mmc.n_pixels_in_row());
    }
    out << static_cast<FileInt>(mmc.n_emissions_);
    out << static_cast<FileInt>(mmc.dr_.detectors());

    for (FileHalf a = 0; a < mmc.dr_.detectors(); ++a) {
      for (FileHalf b = 0; b <= a; ++b) {
        LOR lor(a, b);
        if (mmc.output_triangular) {
          auto pixels = mmc.t_matrix_[t_lor_index(lor)];
          if (pixels) {
            out << a << b;
            // find out count of non-zero pixels
            FileInt count = 0;
            for (auto i = 0; i < mmc.total_n_pixels_in_triangle(); ++i) {
              if (pixels[i])
                count++;
            }
            out << count;

            // write non-zero pixel pairs
            for (FileHalf y = 0; y < mmc.n_pixels_in_row_half(); ++y) {
              for (FileHalf x = 0; x <= y; ++x) {
                FileInt hits = pixels[Super::t_pixel_index(x, y)];
                if (hits) {
                  out << x << y << hits;
                }
              }
            }
          }
        } else {  // output full (may be slow)
                  // find out count of non-zero pixels
          FileInt count = 0;
          for (FileHalf x = 0; x < mmc.n_pixels_in_row(); ++x) {
            for (FileHalf y = 0; y < mmc.n_pixels_in_row(); ++y) {
              if (mmc(lor, x, y))
                count++;
            }
          }
          if (count) {
            out << a << b << count;
            // write out non-zero hits
            for (FileHalf x = 0; x < mmc.n_pixels_in_row(); ++x) {
              for (FileHalf y = 0; y < mmc.n_pixels_in_row(); ++y) {
                FileInt hits = mmc(lor, x, y);
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

  friend ibstream& operator>>(ibstream& in, MatrixLORMajor& mmc) {
    // read header
    FileInt in_is_triangular, in_n_pixels, in_n_emissions, in_n_detectors;
    mmc.read_header(
        in, in_is_triangular, in_n_pixels, in_n_emissions, in_n_detectors);

    if (mmc.n_emissions_ && !in_is_triangular) {
      throw("full matrix cannot be loaded to non-empty matrix");
    }

    // validate incoming parameters
    if (in_n_pixels &&
        in_n_pixels != (in_is_triangular ? mmc.n_pixels_in_row_half() :
                            mmc.n_pixels_in_row())) {
      std::ostringstream msg;
      msg << "incompatible input matrix dimensions " << in_n_pixels << " != "
          << mmc.n_pixels_in_row_half();
      throw(msg.str());
    }
    if (in_n_detectors && in_n_detectors != mmc.dr_.detectors()) {
      throw("incompatible input number of detectors");
    }

    mmc.n_emissions_ += in_n_emissions;

    // load hits
    for (;;) {
      FileHalf a, b;
      in >> a >> b;
      LOR lor(a, b);

      if (in.eof())
        break;

      FileInt count;
      in >> count;

      if (in_is_triangular) {
        auto i_lor = t_lor_index(lor);
        if (i_lor >= mmc.n_lors_) {
          std::ostringstream msg;
          msg << "invalid LOR address (" << a << "," << b << ")";
          throw(msg.str());
        }

        auto pixels = mmc.t_matrix_[i_lor];
        if (!pixels) {
          mmc.t_matrix_[i_lor] = pixels =
                                 new Hit[mmc.total_n_pixels_in_triangle()]();
        }
        // increment hits
        for (auto i = 0; i < count; ++i) {
          FileHalf x, y;
          FileInt hits;
          in >> x >> y >> hits;
          auto i_pixel = Super::t_pixel_index(x, y);
          pixels[i_pixel] += hits;
          mmc.add_hit(i_pixel, hits);
        }
      } else {  // full matrix
                // increment hits
        for (auto i = 0; i < count; ++i) {
          FileHalf x, y;
          FileInt hits;
          in >> x >> y >> hits;

          bool diag;
          int symmetry;
          auto i_pixel = mmc.pixel_index(x, y, diag, symmetry);
          if (i_pixel >= mmc.total_n_pixels_in_triangle())
            throw("invalid pixel address");

          auto i_lor = mmc.lor_index(lor, symmetry);
          if (i_lor >= mmc.n_lors_) {
            std::ostringstream msg;
            msg << "invalid LOR address (" << a << "," << b;
            throw(msg.str());
          }

          auto pixels = mmc.t_matrix_[i_lor];
          if (!pixels) {
            mmc.t_matrix_[i_lor] = pixels =
                                   new Hit[mmc.total_n_pixels_in_triangle()]();
          }
          pixels[i_pixel] += hits;
          mmc.add_hit(i_pixel, hits);
        }
      }
    }

    // we need to divide all values by 8 when loading full matrix
    if (!in_is_triangular) {
      for (auto i_lor = 0; i_lor < mmc.n_lors_; ++i_lor) {
        auto pixels = mmc.t_matrix_[i_lor];
        if (!pixels)
          continue;
        for (auto i_pixel = 0; i_pixel < mmc.total_n_pixels_in_triangle();
             ++i_pixel) {
          pixels[i_pixel] /= 8;
        }
      }
    }

    return in;
  }

  // text output (for validation)
  friend std::ostream& operator<<(std::ostream& out, MatrixLORMajor& mmc) {
    out << "n_pixels_=" << mmc.n_pixels_in_row() << std::endl;
    out << "n_emissions=" << mmc.n_emissions_ << std::endl;
    out << "n_detectors=" << mmc.dr_.detectors() << std::endl;

    for (FileHalf a = 0; a < mmc.dr_.detectors(); ++a) {
      for (FileHalf b = 0; b <= a; ++b) {
        LOR lor(a, b);
        auto pixels = mmc.t_matrix_[t_lor_index(lor)];
        if (pixels) {
          out << "  lor=(" << a << "," << b << ")" << std::endl;
          // find out count of non-zero pixels
          FileInt count = 0;
          for (auto i = 0; i < mmc.total_n_pixels_in_triangle(); ++i) {
            if (pixels[i])
              count++;
          }
          out << "  count=" << count << std::endl;

          // write non-zero pixel pairs
          for (FileHalf y = 0; y < mmc.n_pixels_in_row_half(); ++y) {
            for (FileHalf x = 0; x <= y; ++x) {
              FileInt hits = pixels[Super::t_pixel_index(x, y)];
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

  static void read_header(ibstream& in,
                          FileInt& in_is_triangular,
                          FileInt& in_n_pixels,
                          FileInt& in_n_emissions,
                          FileInt& in_n_detectors) {
    FileInt in_magic;
    in >> in_magic;
    if (in_magic != MAGIC_VERSION_TRIANGULAR &&
        in_magic != MAGIC_VERSION_FULL && in_magic != MAGIC_VERSION_1 &&
        in_magic != MAGIC_VERSION_2) {
      throw("invalid file type format");
    }
    in_is_triangular = (in_magic != MAGIC_VERSION_FULL);

    // load matrix size
    in >> in_n_pixels;

    // load number of emissions
    if (in_magic == MAGIC_VERSION_TRIANGULAR ||
        in_magic == MAGIC_VERSION_FULL || in_magic == MAGIC_VERSION_2) {
      in >> in_n_emissions;
    }

    // load number of detectors
    in_n_detectors = 0;
    if (in_magic == MAGIC_VERSION_TRIANGULAR ||
        in_magic == MAGIC_VERSION_FULL) {
      in >> in_n_detectors;
    }
  }

  template <class FileWriter> void output_lor_bitmap(FileWriter& fw, LOR& lor) {
    fw.template write_header<BitmapPixel>(Super::n_pixels_in_row(),
                                          Super::n_pixels_in_row());
    Hit pixel_max = 0;
    for (auto y = 0; y < Super::n_pixels_in_row(); ++y) {
      for (auto x = 0; x < Super::n_pixels_in_row(); ++x) {
        pixel_max = std::max(pixel_max, (*this)(lor, x, y));
      }
    }
    auto gain = static_cast<double>(std::numeric_limits<BitmapPixel>::max()) /
                pixel_max;
    for (int y = this->n_pixels_in_row() - 1; y >= 0; --y) {
      BitmapPixel row[Super::n_pixels_in_row()];
      for (auto x = 0; x < Super::n_pixels_in_row(); ++x) {
        row[x] =
            std::numeric_limits<BitmapPixel>::max() - gain * (*this)(lor, x, y);
      }
      fw.write_row(row);
    }
  }

  bool output_triangular;

 private:
  static size_t t_lor_index(LOR& lor) {
    if (lor.first < lor.second) {
      std::swap(lor.first, lor.second);
    }
    return lor.first * (lor.first + 1) / 2 + lor.second;
  }

  /// Computes LOR index based on given symmetry (1 out 8)
  /// @param lor      detector number pair
  /// @param symmetry number (0..7)
  size_t lor_index(LOR lor, int symmetry) const {
    if (symmetry & 1) {
      lor.first = (n_2_detectors_ - lor.first) % dr_.detectors();
      lor.second = (n_2_detectors_ - lor.second) % dr_.detectors();
    }
    if (symmetry & 2) {
      lor.first = (n_1_detectors_2_ - lor.first) % dr_.detectors();
      lor.second = (n_1_detectors_2_ - lor.second) % dr_.detectors();
    }
    if (symmetry & 4) {
      lor.first = (n_1_detectors_4_ - lor.first) % dr_.detectors();
      lor.second = (n_1_detectors_4_ - lor.second) % dr_.detectors();
    }
    return t_lor_index(lor);
  }

  DetectorRing<F>& dr_;
  Matrix t_matrix_;
  Hit n_emissions_;
  size_t n_2_detectors_;
  size_t n_1_detectors_2_;
  size_t n_1_detectors_4_;
  size_t n_lors_;
  F s_pixel_;
};
