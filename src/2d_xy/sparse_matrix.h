/// \mainpage
///
/// \section Sparse system matrix binary file format
///
///   - \c uint32_t \b magic =
///     -# \c PETp  triangular
///     -# \c PETP  full
///     -# \c TOFp  TOF triangular
///     -# \c TOFP  TOF full
///   - \c uint32_t \b n_pixels_    half size for \c PETp,
///                                 full size for \c PETP
///   - \c uint32_t \b n_emissions  per pixel
///   - \c uint32_t \b n_detectors  regardless of magic
///   - \c while (!eof)
///     - \c uint16_t \b lor_a, \b lor_b  pair
///     - \c uint32_t \b pixel_pair_count
///     - \c for(.. count ..)
///       - \c uint32_t \b position             only for TOF type
///       - \c uint16_t \b pixel_x, \b pixel_y  half pixels for \c PETp,
///                                             pixels for \c PETP
///       - \c uint32_t \b pixel_hits
///
/// \b Note: TOF position has no particular meaning without quantisation
/// definition. However this is not needed for reconstruction.

#pragma once

#include "util/bstream.h"

#include <vector>

#define fourcc(str, sym) static const uint32_t sym = (*(int*)(str))

// binary serialization                        //  pixels detectors triangular
fourcc("PETt", MAGIC_VERSION_1);               //                       X
fourcc("PETs", MAGIC_VERSION_2);               //     X                 X
fourcc("PETp", MAGIC_VERSION_TRIANGULAR);      //     X        X        X
fourcc("PETP", MAGIC_VERSION_FULL);            //     X        X
fourcc("TOFp", MAGIC_VERSION_TOF_TRIANGULAR);  //     X        X        X
fourcc("TOFP", MAGIC_VERSION_TOF_FULL);        //     X        X

template <typename LORType,
          typename PositionType,
          typename PixelType,
          typename HitType = int>
struct SparseElement {
  SparseElement(LORType&& lor,
                PositionType&& position,
                PixelType&& pixel,
                HitType&& hits)
      : lor(lor), position(position), pixel(pixel), hits(hits) {}
  SparseElement(const LORType& lor,
                const PositionType& position,
                const PixelType& pixel,
                const HitType& hits)
      : lor(lor), position(position), pixel(pixel), hits(hits) {}
  SparseElement() = default;

  LORType lor;
  PositionType position;
  PixelType pixel;
  HitType hits;
};

template <typename PixelType,
          typename LORType,
          typename SType = int,
          typename HitType = int>
class SparseMatrix
    : public std::vector<SparseElement<LORType, SType, PixelType, HitType>> {
  typedef std::vector<SparseElement<LORType, SType, PixelType, HitType>> Super;

 public:
  typedef PixelType Pixel;
  typedef LORType LOR;
  typedef SType S;
  typedef typename std::make_signed<S>::type SS;
  typedef HitType Hit;
  typedef SparseElement<LOR, SType, Pixel, Hit> Element;

  // file representation types, size independent
  typedef uint8_t BitmapPixel;
  typedef uint32_t FileInt;
  typedef uint16_t FileHalf;

  SparseMatrix(S n_pixels_in_row,
               S n_detectors,
               S n_emissions,
               S n_tof_positions = 1,
               bool triangular = true)
      : n_pixels_in_row_(n_pixels_in_row),
        n_pixels_in_row_half_(n_pixels_in_row / 2),
        n_detectors_(n_detectors),
        n_2_detectors_(2 * n_detectors),
        n_emissions_(n_emissions),
        n_lors_(LOR::end_for_detectors(n_detectors).index()),
        triangular_(triangular),
        n_tof_positions_(n_tof_positions) {
    S n_detectors_2 = n_detectors / 2;
    S n_detectors_4 = n_detectors / 4;
    n_1_detectors_2_ = n_detectors + n_detectors_2;
    n_1_detectors_4_ = n_detectors + n_detectors_4;
  }

  S n_pixels_in_row() const { return n_pixels_in_row_; }
  S n_pixels_in_row_half() const { return n_pixels_in_row_half_; }
  S n_detectors() const { return n_detectors_; }
  S n_emissions() const { return n_emissions_; }
  S n_tof_positions() const { return n_tof_positions_; }
  bool triangular() const { return triangular_; }

  SparseMatrix(ibstream& in) {
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

  SparseMatrix& operator<<(const SparseMatrix& other) {

    if (n_pixels_in_row_ != other.n_pixels_in_row_ ||
        n_detectors_ != other.n_detectors_ ||
        triangular_ != other.triangular_ ||
        n_tof_positions_ != other.n_tof_positions_) {
      throw("cannot join two incompatible sparse matrices");
    }

    n_emissions_ += other.n_emissions_;

    this->reserve(this->size() + other.size());
    this->insert(this->end(), other.begin(), other.end());

    auto& first = this->front();
    sort_by_pixel_n_lor();
    for (auto& e : *this) {
      if (first.lor != e.lor || first.position != e.position ||
          first.pixel != e.pixel) {
        first = e;
      } else {
        first.hits += e.hits;
        e.hits = 0;
      }
    }
    sort_by_pixel_n_lor_leaving_empty();
    auto first_empty = std::lower_bound(Super::begin(),
                                        Super::end(),
                                        Element(LOR(), 0, Pixel(), 0),
                                        SortByPixelNPositionNLORLeavingEmpty());
    this->resize(first_empty - this->begin());
    return *this;
  }

  friend obstream& operator<<(obstream& out, SparseMatrix& sm) {
    auto tof = (sm.n_tof_positions_ > 1);
    if (sm.triangular_) {
      out << (tof ? MAGIC_VERSION_TOF_TRIANGULAR : MAGIC_VERSION_TRIANGULAR);
      out << static_cast<FileInt>(sm.n_pixels_in_row_ / 2);
    } else {
      out << (tof ? MAGIC_VERSION_TOF_FULL : MAGIC_VERSION_FULL);
      out << static_cast<FileInt>(sm.n_pixels_in_row_);
    }
    out << static_cast<FileInt>(sm.n_emissions_);
    out << static_cast<FileInt>(sm.n_detectors_);
    if (tof) {
      out << static_cast<FileInt>(sm.n_tof_positions_);
    }

    sm.sort_by_lor();

    LOR current_lor = LOR::end_for_detectors(sm.n_detectors_);

    for (auto it = sm.begin(); it != sm.end(); ++it) {
      auto lor = it->lor;
      auto pixel = it->pixel;
      auto hits = it->hits;

      if (lor != current_lor) {
        current_lor = lor;
        // write down LOR
        out << static_cast<FileHalf>(current_lor.first)
            << static_cast<FileHalf>(current_lor.second);
        // find out count of current LOR elements
        FileInt count = 0;
        for (auto cit = it; cit != sm.end(); ++cit, ++count) {
          if (cit->lor != current_lor)
            break;
        }
        out << count;
      }
      if (tof) {
        out << it->position << static_cast<FileHalf>(pixel.x)
            << static_cast<FileHalf>(pixel.y) << hits;
      } else {
        out << static_cast<FileHalf>(pixel.x) << static_cast<FileHalf>(pixel.y)
            << hits;
      }
    }

    return out;
  }

  // text output (for validation)
  friend std::ostream& operator<<(std::ostream& out, SparseMatrix& sm) {
    out << "pixels in row: " << sm.n_pixels_in_row_ << std::endl;
    out << "    emissions: " << sm.n_emissions_ << std::endl;
    out << "    detectors: " << sm.n_detectors_ << std::endl;

    for (auto it = sm.begin(); it != sm.end(); ++it) {
      if (it->lor.first == it->lor.second) {
        std::ostringstream msg;
        msg << __PRETTY_FUNCTION__ << " invalid LOR (" << it->lor.first << ", "
            << it->lor.second << ")";
        throw(msg.str());
      }
      out << " lor: (" << it->lor.first << ", " << it->lor.second << ")"
          << " position: " << it->position << " pixel: (" << it->pixel.x << ","
          << it->pixel.y << ")"
          << " hits: " << it->hits << std::endl;
    }

    return out;
  }

  template <class FileWriter>
  void output_bitmap(FileWriter& fw, LOR lor = LOR(), S position = -1) {
    S* pixels = new S[n_pixels_in_row_ * n_pixels_in_row_]();
    if (lor.first != lor.second) {
      sort_by_lor();
      for (auto it = std::lower_bound(Super::begin(),
                                      Super::end(),
                                      Element(lor, position, Pixel(), 0),
                                      SortByLORNPosition());
           it->lor == lor && (position < 0 || position == it->position);
           ++it) {
        auto x = it->pixel.x;
        auto y = it->pixel.y;
        if (triangular_) {
          x += n_pixels_in_row_half_;
          y += n_pixels_in_row_half_;
        }
        pixels[n_pixels_in_row_ * y + x] += it->hits;
      }
    } else {
      if (triangular_) {
        for (auto& e : *this) {
          for (auto symmetry = 0; symmetry < 8; ++symmetry) {
            auto pixel = symmetric_pixel(e.pixel, symmetry);
            pixels[n_pixels_in_row_ * pixel.y + pixel.x] += e.hits;
          }
        }
      } else {
        for (auto& e : *this) {
          pixels[n_pixels_in_row_ * e.pixel.y + e.pixel.x] += e.hits;
        }
      }
    }
    fw.template write_header<BitmapPixel>(n_pixels_in_row_, n_pixels_in_row_);
    Hit pixel_max = 0;
    for (auto p = 0; p < n_pixels_in_row_ * n_pixels_in_row_; ++p) {
      pixel_max = std::max(pixel_max, pixels[p]);
    }
    auto gain =
        pixel_max > 0
            ? static_cast<double>(std::numeric_limits<BitmapPixel>::max()) /
                  pixel_max
            : 0.;
    for (SS y = n_pixels_in_row_ - 1; y >= 0; --y) {
      BitmapPixel row[n_pixels_in_row_];
      for (auto x = 0; x < n_pixels_in_row_; ++x) {
        row[x] = std::numeric_limits<BitmapPixel>::max() -
                 gain * pixels[n_pixels_in_row_ * y + x];
      }
      fw.write_row(row);
    }
  }

  void sort_by_lor() {
    if (n_tof_positions_ > 1) {
      std::sort(Super::begin(), Super::end(), SortByLORNPosition());
    } else {
      std::sort(Super::begin(), Super::end(), SortByLOR());
    }
  }
  void sort_by_pixel() {
    std::sort(Super::begin(), Super::end(), SortByPixel());
  }
  void sort_by_lor_n_pixel() {
    std::sort(Super::begin(), Super::end(), SortByLORNPositionNPixel());
  }
  void sort_by_pixel_n_lor() {
    std::sort(Super::begin(), Super::end(), SortByPixelNPositionNLOR());
  }
  void sort_by_pixel_n_lor_leaving_empty() {
    std::sort(
        Super::begin(), Super::end(), SortByPixelNPositionNLORLeavingEmpty());
  }

  SparseMatrix to_full() {
    if (!triangular_) {
      return *this;
    }
    SparseMatrix full(
        n_pixels_in_row_, n_detectors_, n_emissions_, n_tof_positions_, false);
    full.reserve(this->size() * 8);
    for (auto& e : *this) {
      for (auto symmetry = 0; symmetry < 8; ++symmetry) {
        auto pixel = e.pixel;
        auto hits = e.hits;
// FIXME: this is not valid solution below, but converting to
// full matrix we likely get two entries for same pixel, but this does
// not hurt reconstruction though.
#if 0
        // check if we are at diagonal
        if (pixel.x == pixel.y) {
          // avoid writing diagonals twice
          if (symmetry & 4)
            continue;
          // pixels at diagonal get only half of entries
          hits *= 2;
        }
#endif
        auto lor = symmetric_lor(e.lor, symmetry);
        auto position = e.position;
        // if LOR is swapped, then position should be too
        if (lor.first < lor.second) {
          std::swap(lor.first, lor.second);
          // position should be adjusted here so it always goes from
          // higher detector index to lower
          position = n_tof_positions_ - 1 - position;
        }
        full.push_back(
            Element(lor, position, symmetric_pixel(pixel, symmetry), hits));
      }
    }
    return full;
  }

 private:
  S n_pixels_in_row_;
  S n_pixels_in_row_half_;
  S n_detectors_;
  S n_2_detectors_;
  S n_1_detectors_2_;
  S n_1_detectors_4_;
  S n_emissions_;
  S n_lors_;
  bool triangular_;
  S n_tof_positions_;

  struct SortByPixel {
    bool operator()(const Element& a, const Element& b) const {
      return a.pixel < b.pixel;
    }
  };

  struct SortByLOR {
    bool operator()(const Element& a, const Element& b) const {
      return a.lor < b.lor;
    }
  };

#define SparseMatrixCompareField(a, b, field) \
  if (a.field < b.field) {                    \
    return true;                              \
  }                                           \
  if (a.field > b.field) {                    \
    return false;                             \
  }

  struct SortByPixelNPositionNLOR {
    bool operator()(const Element& a, const Element& b) const {
      SparseMatrixCompareField(a, b, pixel.y);
      SparseMatrixCompareField(a, b, pixel.x);
      SparseMatrixCompareField(a, b, position);
      SparseMatrixCompareField(a, b, lor);
      return false;
    }
  };

  struct SortByPixelNPositionNLORLeavingEmpty {
    bool operator()(const Element& a, const Element& b) const {
      if (a.hits && !b.hits)
        return true;
      if (!a.hits && b.hits)
        return false;
      SparseMatrixCompareField(a, b, pixel.y);
      SparseMatrixCompareField(a, b, pixel.x);
      SparseMatrixCompareField(a, b, position);
      SparseMatrixCompareField(a, b, lor);
      return false;
    }
  };

  struct SortByLORNPosition {
    bool operator()(const Element& a, const Element& b) const {
      SparseMatrixCompareField(a, b, lor);
      SparseMatrixCompareField(a, b, position);
      return false;
    }
  };

  struct SortByLORNPositionNPixel {
    bool operator()(const Element& a, const Element& b) const {
      SparseMatrixCompareField(a, b, lor);
      SparseMatrixCompareField(a, b, position);
      SparseMatrixCompareField(a, b, pixel.y);
      SparseMatrixCompareField(a, b, pixel.x);
      return false;
    }
  };

  /// Computes LOR based on given symmetry (1 out 8)
  /// @param lor      base lor for symmetry
  /// @param symmetry number (0..7)
  LOR symmetric_lor(LOR lor, S symmetry) const {
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
    return lor;
  }

  Pixel symmetric_pixel(Pixel p, S symmetry) const {
    if (symmetry & 2) {
      p.x = -p.x - 1;
    }
    if (symmetry & 1) {
      p.y = -p.y - 1;
    }
    // triangulate
    if (symmetry & 4) {
      std::swap(p.x, p.y);
    }
    p.x += n_pixels_in_row_half();
    p.y += n_pixels_in_row_half();
    return p;
  }
};
