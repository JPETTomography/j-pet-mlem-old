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

#include "bstream.h"

template <typename LORType,
          typename PositionType,
          typename PixelType,
          typename HitType = int>
struct SparseElement {
  SparseElement(LORType && lor_a,
                PositionType && position_a,
                PixelType && pixel_a,
                HitType && hits_a)
      : lor(lor_a),
        position(position_a),
        pixel(pixel_a),
        hits(hits_a) {
  }
  SparseElement(const LORType& lor_a,
                const PositionType& position_a,
                const PixelType& pixel_a,
                const HitType& hits_a)
      : lor(lor_a),
        position(position_a),
        pixel(pixel_a),
        hits(hits_a) {
  }
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
  typedef std::vector<
      SparseElement<LORType, SType, PixelType, HitType>> Super;

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

  SparseMatrix(
      S n_pixels_in_row, S n_detectors, S n_emissions, bool triangular = false)
      : n_pixels_in_row_(n_pixels_in_row),
        n_pixels_in_row_half_(n_pixels_in_row / 2),
        n_detectors_(n_detectors),
        n_2_detectors_(2 * n_detectors),
        n_emissions_(n_emissions),
        n_lors_(LOR::end_for_detectors(n_detectors).index()),
        triangular_(triangular) {
    S n_detectors_2 = n_detectors / 2;
    S n_detectors_4 = n_detectors / 4;
    n_1_detectors_2_ = n_detectors + n_detectors_2;
    n_1_detectors_4_ = n_detectors + n_detectors_4;
  }

  S n_pixels_in_row() const { return n_pixels_in_row_; }
  S n_pixels_in_row_half() const { return n_pixels_in_row_half_; }
  S n_detectors() const { return n_detectors_; }
  S n_emissions() const { return n_emissions_; }
  bool triangular() const { return triangular_; }

  SparseMatrix(ibstream& in) {
    FileInt in_magic;
    in >> in_magic;
    if (in_magic != MAGIC_VERSION_TRIANGULAR &&
        in_magic != MAGIC_VERSION_FULL && in_magic != MAGIC_VERSION_1 &&
        in_magic != MAGIC_VERSION_2) {
      throw("invalid file type format");
    }

    FileInt in_is_triangular = (in_magic != MAGIC_VERSION_FULL);

    // load matrix size
    FileInt in_n_pixels_in_row;
    in >> in_n_pixels_in_row;
    if (in_magic != MAGIC_VERSION_FULL)
      in_n_pixels_in_row *= 2;

    // load number of emissions
    FileInt in_n_emissions = 0;
    if (in_magic == MAGIC_VERSION_TRIANGULAR ||
        in_magic == MAGIC_VERSION_FULL || in_magic == MAGIC_VERSION_2) {
      in >> in_n_emissions;
    }

    // load number of detectors
    FileInt in_n_detectors = 0;
    if (in_magic == MAGIC_VERSION_TRIANGULAR ||
        in_magic == MAGIC_VERSION_FULL) {
      in >> in_n_detectors;
    }

    triangular_ = in_is_triangular;
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
        FileInt hits;
        in >> x >> y >> hits;

        this->push_back(Element(lor, 0, Pixel(x, y), hits));
      }
    }
  }

  friend obstream& operator<<(obstream& out, SparseMatrix& sm) {
    if (sm.triangular_) {
      out << MAGIC_VERSION_TRIANGULAR;
      out << static_cast<FileInt>(sm.n_pixels_in_row_ / 2);
    } else {
      out << MAGIC_VERSION_FULL;
      out << static_cast<FileInt>(sm.n_pixels_in_row_);
    }
    out << static_cast<FileInt>(sm.n_emissions_);
    out << static_cast<FileInt>(sm.n_detectors_);

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
      out << static_cast<FileHalf>(pixel.x) << static_cast<FileHalf>(pixel.y)
          << hits;
    }

    return out;
  }

  // text output (for validation)
  friend std::ostream& operator<<(std::ostream& out, SparseMatrix& sm) {
    out << "pixels in row: " << sm.n_pixels_in_row_ << std::endl;
    out << "    emissions: " << sm.n_emissions_ << std::endl;
    out << "    detectors: " << sm.n_detectors_ << std::endl;

    for (auto it = sm.begin(); it != sm.end(); ++it) {
      out << " lor: (" << it->lor.first << ", " << it->lor.second << ")"
          << " pixel: (" << it->pixel.x << "," << it->pixel.y << ")"
          << " hits: " << it->hits << std::endl;
    }

    return out;
  }

  template <class FileWriter>
  void output_lor_bitmap(FileWriter& fw __attribute__((unused)),
                         LOR& lor __attribute__((unused))) {
#if 0
    // FIXME: implement me!
    fw.template write_header<BitmapPixel>(n_pixels_in_row_, n_pixels_in_row_);
    Hit pixel_max = 0;
    for (auto y = 0; y < n_pixels_in_row_; ++y) {
      for (auto x = 0; x < n_pixels_in_row_; ++x) {
        pixel_max = std::max(pixel_max, (*this)(lor, x, y));
      }
    }
    auto gain = static_cast<double>(std::numeric_limits<BitmapPixel>::max()) /
                pixel_max;
    for (SS y = n_pixels_in_row_ - 1; y >= 0; --y) {
      BitmapPixel row[n_pixels_in_row_];
      for (auto x = 0; x < n_pixels_in_row_; ++x) {
        row[x] =
            std::numeric_limits<BitmapPixel>::max() - gain * (*this)(lor, x, y);
      }
      fw.write_row(row);
    }
#endif
  }

  void sort_by_lor() { std::sort(Super::begin(), Super::end(), SortByLOR()); }
  void sort_by_pixel() {
    std::sort(Super::begin(), Super::end(), SortByPixel());
  }

  SparseMatrix to_full() {
    if (!triangular_) {
      return *this;
    }
    SparseMatrix full(n_pixels_in_row_, n_detectors_, n_emissions_);
    full.reserve(this->size() * 8);
    for (auto it = this->begin(); it != this->end(); ++it) {
      for (auto symmetry = 0; symmetry < 8; ++symmetry) {
        auto pixel = it->pixel;
        // avoid writing diagonals twice
        if (symmetry & 4 && pixel.x == pixel.y)
          continue;
        full.push_back(Element(symmetric_lor(it->lor, symmetry),
                               0,
                               symmetric_pixel(pixel, symmetry),
                               it->hits));
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

#define fourcc(a, b, c, d) (((d) << 24) | ((c) << 16) | ((b) << 8) | (a))

 public:
  // binary serialization          // n_pixels_  n_detectors  triagular
  static const FileInt MAGIC_VERSION_1 =
      fourcc('P', 'E', 'T', 't');  //                           X
  static const FileInt MAGIC_VERSION_2 =
      fourcc('P', 'E', 'T', 's');  //     X                     X
  static const FileInt MAGIC_VERSION_TRIANGULAR =
      fourcc('P', 'E', 'T', 'p');  //     X          X          X
  static const FileInt MAGIC_VERSION_FULL =
      fourcc('P', 'E', 'T', 'P');  //     X          X
};
