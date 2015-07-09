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
template <typename LORType,
          typename PositionType,
          typename PixelType,
          typename HitType>
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

/// Made for efficient storage of large PET system matrix.
/// \image html detector_ring2.pdf.png
/// \see \ref sparse_format
template <typename PixelType,
          typename LORType,
          typename SType,
          typename HitType>
class SparseMatrix
    : public std::vector<SparseElement<LORType, SType, PixelType, HitType>> {
  using Base = std::vector<SparseElement<LORType, SType, PixelType, HitType>>;

 public:
  using Pixel = PixelType;
  using LOR = LORType;
  using S = SType;
  using SS = typename std::make_signed<S>::type;
  using Hit = HitType;
  using Element = SparseElement<LOR, SType, Pixel, Hit>;

  // file representation types, size independent
  using BitmapPixel = uint8_t;
  using FileInt = uint32_t;
  using FileHalf = uint16_t;

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

  SparseMatrix(S n_pixels_in_row,
               S n_detectors,
               Hit n_emissions,
               S n_tof_positions = 1,
               bool triangular = true)
      : n_pixels_in_row_(n_pixels_in_row),
        n_pixels_in_row_half_(n_pixels_in_row / 2),
        n_detectors_(n_detectors),
        n_emissions_(n_emissions),
        n_lors_(LOR::end_for_detectors(n_detectors).index()),
        triangular_(triangular),
        n_tof_positions_(n_tof_positions) {}

  S n_pixels_in_row() const { return n_pixels_in_row_; }
  S n_pixels_in_row_half() const { return n_pixels_in_row_half_; }
  S n_detectors() const { return n_detectors_; }
  Hit n_emissions() const { return n_emissions_; }
  S n_tof_positions() const { return n_tof_positions_; }
  bool triangular() const { return triangular_; }

  SparseMatrix(util::ibstream& in) {
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
  };

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
    return *this;
  }

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

  void sort_by_lor() {
    if (n_tof_positions_ > 1) {
      std::sort(Base::begin(), Base::end(), SortByLORNPosition());
    } else {
      std::sort(Base::begin(), Base::end(), SortByLOR());
    }
  }
  void sort_by_pixel() { std::sort(Base::begin(), Base::end(), SortByPixel()); }
  void sort_by_lor_n_pixel() {
    std::sort(Base::begin(), Base::end(), SortByLORNPositionNPixel());
  }
  void sort_by_pixel_n_lor() {
    std::sort(Base::begin(), Base::end(), SortByPixelNPositionNLOR());
  }
  void sort_by_pixel_n_lor_leaving_empty() {
    std::sort(
        Base::begin(), Base::end(), SortByPixelNPositionNLORLeavingEmpty());
  }

  SparseMatrix to_full(
      PET2D::Barrel::SymmetryDescriptor<S>* symmetry_descriptor) {
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
        auto lor = LOR(
            symmetry_descriptor->symmetric_detector(e.lor.first, symmetry),
            symmetry_descriptor->symmetric_detector(e.lor.second, symmetry));
        auto position = e.position;
        // if LOR is swapped, then position should be too
        if (lor.first < lor.second) {
          std::swap(lor.first, lor.second);
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
