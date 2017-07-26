#pragma once

#include "util/bstream.h"
#include <iostream>
#include <vector>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <limits>
#include <algorithm>

#include "2d/geometry/pixel_grid.h"
#include "2d/barrel/symmetry_descriptor.h"

namespace PET2D {
namespace Barrel {

/// Single element of sparse PET system matrix
template <typename LORType, typename PixelType, typename HitType>
struct SparseElement {
  using LOR = LORType;
  using Position = typename LOR::S;
  using Pixel = PixelType;
  using Hit = HitType;

  SparseElement(LOR&& lor, Position&& position, Pixel&& pixel, Hit&& hits)
      : lor(lor), position(position), pixel(pixel), hits(hits) {}

  SparseElement(const LOR& lor,
                const Position& position,
                const Pixel& pixel,
                const Hit& hits)
      : lor(lor), position(position), pixel(pixel), hits(hits) {}

  SparseElement() = default;

  LOR lor;
  Position position;
  Pixel pixel;
  Hit hits;
};

/// \page sparse_format Sparse system matrix binary file format
///
/// \brief Describes binary format used to keep optimal representation for
/// system matrix.
///
///   - \c uint32_t \b magic in
///     -# \c PETp  triangular
///     -# \c PETP  full
///     -# \c TOFp  TOF triangular
///     -# \c TOFP  TOF full
///   - \c uint32_t \b n_pixels     half size for \c PETp,
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
///
/// \image html detector_ring2.pdf.png

/// Sparse 2D barrel PET system matrix
////
/// Made for efficient storage of large PET system matrix.
/// \image html detector_ring2.pdf.png
/// \see \ref sparse_format
/// \see PET2D::Barrel::Geometry
template <typename PixelType, typename LORType, typename HitType>
class SparseMatrix
    : public std::vector<SparseElement<LORType, PixelType, HitType>> {
  using Base = std::vector<SparseElement<LORType, PixelType, HitType>>;

 public:
  using Pixel = PixelType;
  using LOR = LORType;
  using S = typename std::common_type<typename Pixel::S, typename LOR::S>::type;
  using SS = typename std::make_signed<S>::type;
  using Hit = HitType;
  using Element = SparseElement<LOR, Pixel, Hit>;

  // file representation types, size independent
  using BitmapPixel = uint8_t;  ///< output bitmap pixel value type
  using FileInt = uint32_t;     ///< file integer
  using FileHalf = uint16_t;    ///< file half-size integer

  /// First file word magic value representing type of matrix
  enum Magic : FileInt {
    // binary serialization                //  pixels detectors triangular
    // clang-format off
    VERSION_1               = "PETt"_4cc,  //                       X
    VERSION_2               = "PETs"_4cc,  //     X                 X
    VERSION_TRIANGULAR      = "PETp"_4cc,  //     X        X        X
    VERSION_FULL            = "PETP"_4cc,  //     X        X
    VERSION_TOF_TRIANGULAR  = "TOFp"_4cc,  //     X        X        X
    VERSION_TOF_FULL        = "TOFP"_4cc,  //     X        X
    // clang-format on
  };

  enum Sort {
    UNSORTED = 0,
    BY_LOR,
    BY_PIXEL,
    BY_LOR_N_PIXEL,
    BY_PIXEL_N_LOR,
    BY_PIXEL_N_LOR_LEAVING_EMPTY,
  };

  /// Construct new sparse matrix with given parameters.
  SparseMatrix(S n_pixels_in_row,      ///< number of pixels in row
               S n_detectors,          ///< total number of detectors
               S n_tof_positions = 1,  ///< number of TOF positions (1-no TOF)
               Hit n_emissions = 0,    ///< number of future emissions
               bool triangular = true  ///< is this matrix triangular
               )
      : n_pixels_in_row_(n_pixels_in_row),
        n_pixels_in_row_half_(n_pixels_in_row / 2),
        n_detectors_(n_detectors),
        n_emissions_(n_emissions),
        n_lors_(LOR::end_for_detectors(n_detectors).index()),
        triangular_(triangular),
        n_tof_positions_(n_tof_positions),
        sorted_(UNSORTED) {}

  /// Returns matrix current sort method.
  Sort sorted() const { return sorted_; }

  /// Returns number of pixels in row.
  S n_pixels_in_row() const { return n_pixels_in_row_; }

  /// Returns half number of pixels in row.
  ////
  /// For triangular matrix, this is actual number of pixels covered by it,
  /// since we just keep 1/8 of full matrix.
  S n_pixels_in_row_half() const { return n_pixels_in_row_half_; }

  /// Returns number of detectors.
  S n_detectors() const { return n_detectors_; }

  /// Returns number of emissions described by this matrix.
  Hit n_emissions() const { return n_emissions_; }

  /// Returns number of TOF positions (1-means there is no TOF imformation).
  S n_tof_positions() const { return n_tof_positions_; }

  /// True if matrix is triangular (keeps only 1/8 of full matrix information).
  bool triangular() const { return triangular_; }

  /// Increment number of emissions by given number.
  ////
  /// This is intended when some external Monte-Carlo appends simulation results
  /// manually via \c push_back or \c emplace_back.
  void increment_n_emissions(Hit increment) { n_emissions_ += increment; }

  SparseMatrix(util::ibstream& in) : sorted_(BY_LOR_N_PIXEL) {
    FileInt in_magic;
    in >> in_magic;
    if (in_magic != Magic::VERSION_TRIANGULAR &&
        in_magic != Magic::VERSION_FULL &&
        in_magic != Magic::VERSION_TOF_TRIANGULAR &&
        in_magic != Magic::VERSION_TOF_FULL && in_magic != Magic::VERSION_1 &&
        in_magic != Magic::VERSION_2) {
      throw("invalid file type format");
    }

    bool in_is_triangular = (in_magic != Magic::VERSION_FULL &&
                             in_magic != Magic::VERSION_TOF_FULL);
    bool in_is_tof = (in_magic == Magic::VERSION_TOF_TRIANGULAR ||
                      in_magic == Magic::VERSION_TOF_FULL);

    FileInt in_n_pixels_in_row;
    in >> in_n_pixels_in_row;
    if (in_is_triangular)
      in_n_pixels_in_row *= 2;
    else
      sorted_ = UNSORTED;

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

        this->emplace_back(lor, position, Pixel(x, y), hits);
      }
    }
  }

  /// Merge duplicate hit entries pointing to same LOR-pixel-position.
  /// Matrix has to be sorted first.
  void merge_duplicates() {
    for (auto it_elem = this->begin(); it_elem != this->end(); ++it_elem) {
      auto next = it_elem + 1;
      if (next != this->end()) {
        auto first = *it_elem;
        auto second = *next;
        if (first.lor == second.lor && first.pixel == second.pixel &&
            first.position == second.position) {

          it_elem->hits += second.hits;
          next->hits = 0;
        }
      }
    }

    this->erase(
        std::remove_if(this->begin(),
                       this->end(),
                       [](const Element& a) -> bool { return a.hits == 0; }),
        this->end());
  }

  /// Append other sparse matrix to this matrix.
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
    auto first_empty = std::lower_bound(Base::begin(),
                                        Base::end(),
                                        Element(LOR(), 0, Pixel(), 0),
                                        SortByPixelNPositionNLORLeavingEmpty());
    this->resize(first_empty - this->begin());
    sorted_ = BY_PIXEL_N_LOR;
    return *this;
  }

  /// Output sparse matrix to binary stream.
  friend util::obstream& operator<<(util::obstream& out, SparseMatrix& sm) {
    auto tof = (sm.n_tof_positions_ > 1);
    if (sm.triangular_) {
      out << (tof ? Magic::VERSION_TOF_TRIANGULAR : Magic::VERSION_TRIANGULAR);
      out << static_cast<FileInt>(sm.n_pixels_in_row_ / 2);
    } else {
      out << (tof ? Magic::VERSION_TOF_FULL : Magic::VERSION_FULL);
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
        out << static_cast<FileInt>(it->position)  //
            << static_cast<FileHalf>(pixel.x) << static_cast<FileHalf>(pixel.y)
            << static_cast<FileInt>(hits);
      } else {
        out << static_cast<FileHalf>(pixel.x) << static_cast<FileHalf>(pixel.y)
            << static_cast<FileInt>(hits);
      }
    }

    return out;
  }

  /// Output sparse matrix text stream.
  ////
  /// Used for validation only.
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

  /// Output bitmap of sentitivity or lor (if non empty) to bitmap writer.
  template <class FileWriter>
  void output_bitmap(FileWriter& fw, LOR lor = LOR(), S position = -1) {
    Hit* pixels = new Hit[n_pixels_in_row_ * n_pixels_in_row_]();
    if (lor.first != lor.second) {
      sort_by_lor();
      for (auto it = std::lower_bound(Base::begin(),
                                      Base::end(),
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
    BitmapPixel* row =
        (BitmapPixel*)alloca(n_pixels_in_row_ * sizeof(BitmapPixel));
    for (SS y = n_pixels_in_row_ - 1; y >= 0; --y) {
      for (auto x = 0; x < n_pixels_in_row_; ++x) {
        row[x] = std::numeric_limits<BitmapPixel>::max() -
                 gain * pixels[n_pixels_in_row_ * y + x];
      }
      fw.write_row(row);
    }
  }

  /// Sort matrix entries by LOR index.
  void sort_by_lor() {
    if (sorted_ == BY_LOR_N_PIXEL || sorted_ == BY_LOR)
      return;
    if (n_tof_positions_ > 1) {
      std::sort(Base::begin(), Base::end(), SortByLORNPosition());
    } else {
      std::sort(Base::begin(), Base::end(), SortByLOR());
    }
    sorted_ = BY_LOR;
  }
  /// Sort matrix entries by pixel index.
  void sort_by_pixel() {
    if (sorted_ == BY_PIXEL_N_LOR || sorted_ == BY_PIXEL)
      return;
    std::sort(Base::begin(), Base::end(), SortByPixel());
    sorted_ = BY_PIXEL;
  }
  /// Sort matrix entries by LOR then pixel index.
  void sort_by_lor_n_pixel() {
    if (sorted_ == BY_LOR_N_PIXEL)
      return;
    std::sort(Base::begin(), Base::end(), SortByLORNPositionNPixel());
    sorted_ = BY_LOR_N_PIXEL;
  }
  /// Sort matrix entries by pixel then LOR index.
  void sort_by_pixel_n_lor() {
    if (sorted_ == BY_PIXEL_N_LOR)
      return;
    std::sort(Base::begin(), Base::end(), SortByPixelNPositionNLOR());
    sorted_ = BY_PIXEL_N_LOR;
  }
  /// Sort matrix entries by pixel then LOR index, leaving empty elements back.
  void sort_by_pixel_n_lor_leaving_empty() {
    if (sorted_ == BY_PIXEL_N_LOR_LEAVING_EMPTY)
      return;
    std::sort(
        Base::begin(), Base::end(), SortByPixelNPositionNLORLeavingEmpty());
    sorted_ = BY_PIXEL_N_LOR_LEAVING_EMPTY;
  }

  /// Return full (non-triangular) sparse matrix.
  SparseMatrix to_full(
      const PET2D::Barrel::SymmetryDescriptor<S>& symmetry_descriptor) {
    if (!triangular_) {
      return *this;
    }
    SparseMatrix full(
        n_pixels_in_row_, n_detectors_, n_tof_positions_, n_emissions_, false);
    full.reserve(this->size() * 8);
    for (auto& e : *this) {
      for (auto symmetry = 0; symmetry < 8; ++symmetry) {
        auto pixel = e.pixel;
        auto hits = e.hits;

        // FIXME: The solution below is not valid, but converting to full matrix
        // we likely get two entries for same pixel, but this does not hurt
        // reconstruction though.
        //
        // NOTE: Monte-Carlo implementations ensure that pixels on diagonal get
        // only half of entries, because random emissions points inside diagonal
        // pixels that overflow diagonal get discarded.
        ;
#if HANDLE_DIAGONALS_SPECIALLY
        // check if we are at diagonal
        if (pixel.x == pixel.y) {
          // avoid writing diagonals twice
          if (symmetry & 4)
            continue;
          // pixels at diagonal get only half of entries
          hits *= 2;
        }
#endif
        auto symmetric_lor_first =
            symmetry_descriptor.symmetric_detector(e.lor.first, symmetry);
        auto symmetric_lor_second =
            symmetry_descriptor.symmetric_detector(e.lor.second, symmetry);
        auto lor = LOR(symmetric_lor_first, symmetric_lor_second);
        auto position = e.position;
        // if LOR is swapped, then position should be too
        if (symmetric_lor_first < symmetric_lor_second) {
          // position should be adjusted here so it always goes from
          // higher detector index to lower
          position = n_tof_positions_ - 1 - position;
        }
        full.emplace_back(
            lor, position, symmetric_pixel(pixel, symmetry), hits);
      }
    }
    return full;
  }

 private:
  S n_pixels_in_row_;
  S n_pixels_in_row_half_;
  S n_detectors_;
  Hit n_emissions_;
  int n_lors_;
  bool triangular_;
  S n_tof_positions_;
  Sort sorted_;

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

 public:
  /// Return symmetrix pixel for given symmetry index (symmetry mask).
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
}  // Barrel
}  // PET2D
