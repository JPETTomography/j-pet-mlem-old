#pragma once

#include "bstream.h"
#include "pixel.h"

template <typename LORType, typename SType = int, typename HitType = int>
class SparseMatrix
    : public std::vector<std::tuple<LORType, Pixel<SType>, HitType>> {
  typedef std::vector<std::tuple<LORType, Pixel<SType>, HitType>> Super;

 public:
  typedef ::Pixel<SType> Pixel;
  typedef LORType LOR;
  typedef SType S;
  typedef typename std::make_signed<S>::type SS;
  typedef HitType Hit;
  typedef std::tuple<LOR, Pixel, Hit> Element;

  // file representation types, size independent
  typedef uint8_t BitmapPixel;
  typedef uint32_t FileInt;
  typedef uint16_t FileHalf;

  SparseMatrix(S n_pixels_in_row, S n_detectors, bool triangular = false)
      : n_pixels_in_row_(n_pixels_in_row),
        n_detectors_(n_detectors),
        n_lors_(LOR::end_for_detectors(n_detectors).index()),
        triangular_(triangular) {
  }

  SparseMatrix(ibstream& in) {
    read_header(in, triangular_, n_pixels_in_row_, n_emissions_, n_detectors_);
    n_lors_ = LOR::end_for_detectors(n_detectors_).t_index();

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
      for (auto i = 0; i < count; ++i) {
        FileHalf x, y;
        FileInt hits;
        in >> x >> y >> hits;
        this->push_back(Element(lor, x, y, hits));
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
      Element element = *it;
      LOR& lor = std::get<0>(element);
      Pixel& p = std::get<1>(element);
      FileInt hits = std::get<2>(element);

      if (lor != current_lor) {
        current_lor = lor;
        // write down LOR
        out << static_cast<FileHalf>(current_lor.first)
            << static_cast<FileHalf>(current_lor.second);
        // find out count of current LOR elements
        FileInt count = 0;
        for (auto cit = it; cit != sm.end(); ++cit, ++count) {
          Element element = *cit;
          LOR& lor = std::get<0>(element);
          if (lor != current_lor)
            break;
        }
        out << count;
      }
      out << static_cast<FileHalf>(p.x) << static_cast<FileHalf>(p.y) << hits;
    }

    return out;
  }

  // text output (for validation)
  friend std::ostream& operator<<(std::ostream& out, SparseMatrix& sm) {
    out << "n_pixels_=" << sm.n_pixels_in_row_ << std::endl;
    out << "n_emissions=" << sm.n_emissions_ << std::endl;
    out << "n_detectors=" << sm.n_detectors_ << std::endl;

    LOR current_lor = LOR::end_for_detectors(sm.n_detectors_);

    for (auto it = sm.begin(); it != sm.end(); ++it) {
      Element element = *it;
      LOR& lor = std::get<0>(element);
      Pixel& p = std::get<1>(element);
      Hit hits = std::get<2>(element);

      if (lor != current_lor) {
        current_lor = lor;
        // write down LOR
        out << "  lor="
            << "(" << current_lor.first << "," << current_lor.second << ")"
            << std::endl;
        // find out count of current LOR elements
        S count = 0;
        for (auto cit = it; cit != sm.end(); ++cit, ++count) {
          Element element = *cit;
          LOR& lor = std::get<0>(element);
          if (lor != current_lor)
            break;
        }
        out << "  count=" << count << std::endl;
      }
      out << "    (" << p.x << "," << p.y << ")=" << hits << std::endl;
    }

    return out;
  }

  template <class FileWriter> void output_lor_bitmap(FileWriter& fw, LOR& lor) {
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

 private:
  S n_pixels_in_row_;
  S n_emissions_;
  S n_detectors_;
  S n_lors_;
  bool triangular_;

  Hit operator()(LOR& lor, S x, S y) { throw(__PRETTY_FUNCTION__); }

  struct SortByPixel {
    bool operator()(const Element& a, const Element& b) const {
      auto ax = std::get<1>(a);
      auto bx = std::get<1>(b);
      auto ay = std::get<2>(a);
      auto by = std::get<2>(b);
      return ay < by ? true : ay > by ? false : ax < bx;
    }
  };

  struct SortByLOR {
    bool operator()(const Element& a, const Element& b) const {
      return std::get<0>(a) < std::get<0>(b);
    }
  };

#define fourcc(a, b, c, d) (((d) << 24) | ((c) << 16) | ((b) << 8) | (a))

  // binary serialization                 // n_pixels_  n_detectors  triagular
  static const FileInt MAGIC_VERSION_1 =
      fourcc('P', 'E', 'T', 't');  //                           X
  static const FileInt MAGIC_VERSION_2 =
      fourcc('P', 'E', 'T', 's');  //     X                     X
  static const FileInt MAGIC_VERSION_TRIANGULAR =
      fourcc('P', 'E', 'T', 'p');  //     X          X          X
  static const FileInt MAGIC_VERSION_FULL =
      fourcc('P', 'E', 'T', 'P');  //     X          X

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

};
