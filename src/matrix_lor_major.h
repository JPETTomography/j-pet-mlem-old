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

#include "matrix.h"

#define fourcc(a, b, c, d) (((d) << 24) | ((c) << 16) | ((b) << 8) | (a))

/// Provides storage for 1/8 PET system matrix
template <typename LORType, typename SType = int, typename HitType = int>
class MatrixLORMajor : public Matrix<LORType, SType, HitType> {
  typedef Matrix<LORType, SType, HitType> Super;

 public:
  typedef LORType LOR;
  typedef SType S;
  typedef typename std::make_signed<S>::type SS;
  typedef HitType Hit;
  typedef uint8_t BitmapPixel;
  typedef uint32_t FileInt;
  typedef uint16_t FileHalf;
  typedef Hit* Pixels;
  typedef Pixels* Matrix;

  /// @param n_pixels    number of pixels in each directions
  /// @param n_detectors number of detectors stored in the matrix
  MatrixLORMajor(S n_pixels, S n_detectors)
      : Super(n_pixels, n_detectors),
        output_triangular(true),
        n_emissions_(0),
        n_detectors_(n_detectors),
        n_2_detectors_(2 * n_detectors),
        n_1_detectors_2_(n_detectors + n_detectors / 2),
        n_1_detectors_4_(n_detectors + n_detectors / 4),
        n_lors_(LOR::end_for_detectors(n_detectors).t_index()) {

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

  void add_to_t_matrix(LOR& lor, S i_pixel) {
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

  Hit operator()(LOR lor, S x, S y) const {
    bool diag;
    S symmetry;
    auto i_pixel = this->pixel_index(x, y, diag, symmetry);
    auto pixels = t_matrix_[lor_index(lor, symmetry)];
    return pixels ? pixels[i_pixel] * (diag ? 2 : 1) : 0;
  }

  S non_zero_lors() {
    S non_zero_lors = 0;
    for (auto i = 0; i < n_lors_; ++i) {
      if (t_matrix_[i])
        ++non_zero_lors;
    }
    return non_zero_lors;
  }

  void increase_n_emissions(Hit n_emissions) { n_emissions_ += n_emissions; }

  friend obstream& operator<<(obstream& out, MatrixLORMajor& mmc) {
    if (mmc.output_triangular) {
      out << MAGIC_VERSION_TRIANGULAR;
      out << static_cast<FileInt>(mmc.n_pixels_in_row_half());
    } else {
      out << MAGIC_VERSION_FULL;
      out << static_cast<FileInt>(mmc.n_pixels_in_row());
    }
    out << static_cast<FileInt>(mmc.n_emissions_);
    out << static_cast<FileInt>(mmc.n_detectors_);

    for (FileHalf a = 0; a < mmc.n_detectors_; ++a) {
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
        in_n_pixels != (in_is_triangular ? mmc.n_pixels_in_row_half()
                                         : mmc.n_pixels_in_row())) {
      std::ostringstream msg;
      msg << "incompatible input matrix dimensions " << in_n_pixels << " != "
          << mmc.n_pixels_in_row_half();
      throw(msg.str());
    }
    if (in_n_detectors && in_n_detectors != mmc.n_detectors_) {
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
          S symmetry;
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
    out << "n_detectors=" << mmc.n_detectors_ << std::endl;

    for (FileHalf a = 0; a < mmc.n_detectors_; ++a) {
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
    for (SS y = this->n_pixels_in_row() - 1; y >= 0; --y) {
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
  static S t_lor_index(LOR& lor) {
    if (lor.first < lor.second) {
      std::swap(lor.first, lor.second);
    }
    return lor.first * (lor.first + 1) / 2 + lor.second;
  }

  /// Computes LOR index based on given symmetry (1 out 8)
  /// @param lor      detector number pair
  /// @param symmetry number (0..7)
  S lor_index(LOR lor, S symmetry) const {
    if (symmetry & 1) {
      lor.first = (n_2_detectors_ - lor.first) % n_detectors_;
      lor.second = (n_2_detectors_ - lor.second) % n_detectors_;
    }
    if (symmetry & 2) {
      lor.first = (n_1_detectors_2_ - lor.first) % n_detectors_;
      lor.second = (n_1_detectors_2_ - lor.second) % n_detectors_;
    }
    if (symmetry & 4) {
      lor.first = (n_1_detectors_4_ - lor.first) % n_detectors_;
      lor.second = (n_1_detectors_4_ - lor.second) % n_detectors_;
    }
    return t_lor_index(lor);
  }

  LOR end_;
  Matrix t_matrix_;
  Hit n_emissions_;
  S n_detectors_;
  S n_2_detectors_;
  S n_1_detectors_2_;
  S n_1_detectors_4_;
  S n_lors_;

 public:
  // binary serialization                 // n_pixels_  n_detectors  triagular
  static const FileInt MAGIC_VERSION_1 =
      fourcc('P', 'E', 'T', 't');  //                           X
  static const FileInt MAGIC_VERSION_2 =
      fourcc('P', 'E', 'T', 's');  //     X                     X
  static const FileInt MAGIC_VERSION_TRIANGULAR =
      fourcc('P', 'E', 'T', 'p');  //     X          X          X
  static const FileInt MAGIC_VERSION_FULL =
      fourcc('P', 'E', 'T', 'P');  //     X          X
};
