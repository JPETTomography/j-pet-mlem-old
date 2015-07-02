#pragma once

#include "matrix.h"

namespace PET2D {
namespace Barrel {

/// Pixel-major memory ordered matrix

/// This class represents a system matrix that stores the content in
/// "pixel major" mode. That is for each pixel a list of lors is kept.
///
/// The idea behind it is that the MC simulations are done on per pixel basis.
/// That means that in the same time only few pixels are processed in parallel.
/// On shiva for example the biggest machine has 24 cores.
/// We can afford to alloc memory for lors in an uneficient way providing for
/// quick acces and reduce them after the simulations of this pixel is finished.
///
/// It also means that as different threads are processing different pixels,
/// there is no need for synchronisation in add_pixels.
///
/// The reconstruction is however done using the lor_major matrix.
/// So this class has the possibility to write the matrix down in the
/// triangular lor_major and full lor_major format.

template <typename PixelType,
          typename LORType,
          typename SType,
          typename HitType>
class MatrixPixelMajor : public Matrix<PixelType, LORType, SType, HitType> {
  using Base = Matrix<PixelType, LORType, SType, HitType>;

 public:
  using Pixel = PixelType;
  using LOR = LORType;
  using S = SType;
  using Hit = HitType;
  using SS = typename std::make_signed<S>::type;
  using SparseMatrix = typename Base::SparseMatrix;
  using SparseElement = typename SparseMatrix::Element;

  MatrixPixelMajor(S n_pixels_in_row, S n_detectors, S n_tof_positions = 1)
      : Base(n_pixels_in_row, n_detectors, n_tof_positions),
        n_pixels_in_row_half(n_pixels_in_row / 2),
        n_pixels_(Pixel::end_for_n_pixels_in_row(n_pixels_in_row).index()),
        n_lors(LOR::end_for_detectors(n_detectors).index()),
        size_(0),
        pixel_lor_hits_ptr(new Hit*[n_pixels_]()),
        pixel_lor_hits(n_pixels_),
        pixel_lor_count(n_pixels_),
        index_to_lor(n_lors),
        index_to_pixel(n_pixels_) {
    // store index to LOR mapping
    for (auto lor = this->begin_lor(); lor != this->end_lor(); ++lor) {
      index_to_lor[lor.index()] = lor;
    }
    // store index to pixel mapping
    for (auto pixel = this->begin_pixel(); pixel != this->end_pixel();
         ++pixel) {
      index_to_pixel[pixel.index()] = pixel;
    }
  }

  void hit_lor(const LOR& lor, S position, S i_pixel, Hit hits = 1) {
    if (position >= this->n_tof_positions()) {
      std::ostringstream msg;
      msg << "hit position " << position << " greater than max TOF positions "
          << this->n_tof_positions();
      throw(msg.str());
    }
    if (lor.first == lor.second) {
      std::ostringstream msg;
      msg << __FUNCTION__ << " invalid LOR " << lor.index() << " (" << lor.first
          << ", " << lor.second << ")";
      throw(msg.str());
    }
    if (!pixel_lor_hits_ptr[i_pixel]) {
      pixel_lor_hits_ptr[i_pixel] = new Hit[n_lors * this->n_tof_positions()]();
      // unpack previous values (if any)
      for (auto& e : pixel_lor_hits[i_pixel]) {
        hit_lor(e.lor, e.position, e.pixel.index(), e.hits);
      }
    }

    auto& current_hits =
        pixel_lor_hits_ptr[i_pixel][lor.index() * this->n_tof_positions() +
                                    position];
    if (current_hits == 0) {
      pixel_lor_count[i_pixel]++;
      size_++;
    }
    current_hits += hits;
  }

  ~MatrixPixelMajor() {
    for (int i_pixel = 0; i_pixel < n_pixels_; ++i_pixel) {
      if (pixel_lor_hits_ptr[i_pixel]) {
        delete[] pixel_lor_hits_ptr[i_pixel];
      }
    }
    delete[] pixel_lor_hits_ptr;
  }

  void compact_pixel_index(S i_pixel) {
    if (!pixel_lor_hits_ptr[i_pixel] || i_pixel >= n_pixels_)
      return;

    // ensure we have enough space for the all LORs for that pixel
    pixel_lor_hits[i_pixel].resize(pixel_lor_count[i_pixel]);

    for (int i_lor = 0, lor_count = 0; i_lor < n_lors; ++i_lor) {
      for (S position = 0; position < this->n_tof_positions(); ++position) {
        auto hits =
            pixel_lor_hits_ptr[i_pixel][i_lor * this->n_tof_positions() +
                                        position];
        if (hits > 0) {
          LOR lor = index_to_lor[i_lor];
          if (lor.first == lor.second) {
            std::ostringstream msg;
            msg << __FUNCTION__ << " invalid LOR " << i_lor << " (" << lor.first
                << ", " << lor.second << ") for pixel index " << i_pixel;
            throw(msg.str());
          }
          pixel_lor_hits[i_pixel][lor_count++] = SparseElement(
              index_to_lor[i_lor], position, index_to_pixel[i_pixel], hits);
        }
      }
    }

    delete[] pixel_lor_hits_ptr[i_pixel], pixel_lor_hits_ptr[i_pixel] = NULL;

    std::sort(pixel_lor_hits[i_pixel].begin(),
              pixel_lor_hits[i_pixel].end(),
              SparseElementLORComparator());
  }

  Pixel pixel_at_index(S i_pixel) { return index_to_pixel[i_pixel]; }

  SparseMatrix to_sparse() {
    SparseMatrix sparse(this->n_pixels_in_row(),
                        this->n_detectors(),
                        this->n_emissions(),
                        this->n_tof_positions());
    sparse.reserve(size_);
    for (int i_pixel = 0; i_pixel < n_pixels_; ++i_pixel) {
      for (auto& e : pixel_lor_hits[i_pixel]) {
        sparse.push_back(e);
      }
    }

    return sparse;
  }

  MatrixPixelMajor& operator<<(SparseMatrix& sparse) {
    if (this->n_emissions()) {
      throw("cannot load multiple sparse matrices into pixel-major matrix");
    }

    sparse.sort_by_pixel_n_lor();
    this->add_emissions(sparse.n_emissions());

    int lor_count = 0;
    int i_current_pixel = n_pixels_;

    for (auto& e : sparse) {
      Pixel pixel = e.pixel;
      auto triangular = sparse.triangular();
      if (!triangular) {
        pixel.x -= n_pixels_in_row_half;
        pixel.y -= n_pixels_in_row_half;
        if (pixel.x < 0 || pixel.y < 0 || pixel.y < pixel.x)
          continue;
      }

      // if we are at diagonal in full matrix, we should have there
      // half of entries
      auto hits = e.hits;
      if (!triangular && pixel.x == pixel.y)
        hits /= 2;

      auto i_pixel = pixel.index();
      if (i_current_pixel != i_pixel) {
        if (lor_count) {
          pixel_lor_count[i_current_pixel] += lor_count;
        }
        lor_count = 0;
        i_current_pixel = i_pixel;
      }
      this->hit(i_pixel, hits);
      pixel_lor_hits[i_pixel].emplace_back(e.lor, e.position, pixel, hits);
    }
    if (lor_count) {
      pixel_lor_count[i_current_pixel] += lor_count;
    }
    return *this;
  }

  // for testing purposes
  Hit lor_hits_at_pixel_index(LOR lor, S i_pixel) {
    auto it =
        std::lower_bound(pixel_lor_hits[i_pixel].begin(),
                         pixel_lor_hits[i_pixel].end(),
                         SparseElement(lor, 0, index_to_pixel[i_pixel], 0),
                         SparseElementLORComparator());

    if (it == pixel_lor_hits[i_pixel].end())
      return 0;
    return it->hits;
  }

  S size() const { return size_; }
  int n_lors_at_pixel_index(S i_pixel) const {
    return pixel_lor_count[i_pixel];
  }
  int n_pixels() const { return n_pixels_; }

 private:
  // disable copy contructor
  MatrixPixelMajor(const MatrixPixelMajor& rhs) : Base(0, 0) {
    (void)rhs;  // unused
    throw(__PRETTY_FUNCTION__);
  }

  struct SparseElementLORComparator {
    bool operator()(const SparseElement& a, const SparseElement& b) const {
      return a.lor < b.lor;
    }
  };

  S n_pixels_in_row_half;
  int n_pixels_;
  int n_lors;
  int size_;
  Hit** pixel_lor_hits_ptr;
  std::vector<std::vector<SparseElement>> pixel_lor_hits;
  std::vector<S> pixel_lor_count;
  std::vector<LOR> index_to_lor;
  std::vector<Pixel> index_to_pixel;
};
}  // Barrel
}  // PET2D
