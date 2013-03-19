#pragma once

#include "matrix.h"

/// This class represents a system matrix that stores the content in
/// "pixel major" mode. That is for each pixel a list w lors is kept.
///
/// The idea behind it is that the MC simulations are done on per pixel basis.
/// That means that in the same time only few pixels are processed in parallel.
/// On shiva for example the biggest machine has 24 cores.
/// We can afford to alloc memory for lors in an uneficient way providing for
/// quick acces and reduce them after the simulations of  this pixel is finished.
///
/// It also means that as different threads are processing different pixels,
/// there is no need for synchronisation in add_pixels.
///
/// The reconstruction is however done using the lor_major matrix.
/// So this class has the possibility to write the matrix down in the
/// triangular lor_major and full lor_major format.

template <typename PixelType,
          typename LORType,
          typename SType = int,
          typename HitType = int>
class MatrixPixelMajor : public Matrix<PixelType, LORType, SType, HitType> {
  typedef Matrix<PixelType, LORType, SType, HitType> Super;

 public:
  typedef PixelType Pixel;
  typedef LORType LOR;
  typedef SType S;
  typedef HitType Hit;
  typedef typename std::make_signed<S>::type SS;
  typedef typename Super::SparseMatrix SparseMatrix;
  typedef typename SparseMatrix::Element SparseElement;

  MatrixPixelMajor(S n_pixels_in_row, S n_detectors, S n_tof_positions = 1)
      : Super(n_pixels_in_row, n_detectors, n_tof_positions),
        n_pixels_in_row_half_(n_pixels_in_row / 2),
        n_pixels_(Pixel::end_for_n_pixels_in_row(n_pixels_in_row).index()),
        n_lors_(LOR::end_for_detectors(n_detectors).index()),
        size_(0),
        pixel_lor_hits_ptr_(new S* [n_pixels_]()),
        pixel_lor_hits_(n_pixels_),
        pixel_lor_count_(n_pixels_),
        index_to_lor_(n_lors_),
        index_to_pixel_(n_pixels_) {
    // store index to LOR mapping
    for (auto lor = this->begin_lor(); lor != this->end_lor(); ++lor) {
      index_to_lor_[lor.index()] = lor;
    }
    // store index to pixel mapping
    for (auto pixel = this->begin_pixel(); pixel != this->end_pixel();
         ++pixel) {
      index_to_pixel_[pixel.index()] = pixel;
    }
  }

  void hit_lor(const LOR& lor, S position, S i_pixel, S hits) {
    if (position >= this->n_tof_positions()) {
      std::ostringstream msg;
      msg << "hit position " << position << " greater than max TOF positions"
          << this->n_tof_positions();
      throw(msg.str());
    }
    if (lor.first == lor.second) {
      std::ostringstream msg;
      msg << __PRETTY_FUNCTION__ << "invalid LOR " << lor.index() << " ("
          << lor.first << ", " << lor.second << ")";
      throw(msg.str());
    }
    if (pixel_lor_hits_ptr_[i_pixel] == 0x0) {
      pixel_lor_hits_ptr_[i_pixel] = new S[n_lors_ * this->n_tof_positions()]();
      // unpack previous values (if any)
      for (auto it = pixel_lor_hits_[i_pixel].begin();
           it != pixel_lor_hits_[i_pixel].end();
           ++it) {
        hit_lor(it->lor, it->position, i_pixel, it->hits);
      }
    }

    auto& current_hits = pixel_lor_hits_ptr_[i_pixel][
        lor.index() * this->n_tof_positions() + position];
    if (current_hits == 0) {
      pixel_lor_count_[i_pixel]++;
      size_++;
    }
    current_hits += hits;
  }

  ~MatrixPixelMajor() {
    for (S i_pixel = 0; i_pixel < n_pixels_; ++i_pixel) {
      if (pixel_lor_hits_ptr_[i_pixel]) {
        delete[] pixel_lor_hits_ptr_[i_pixel];
      }
    }
    delete[] pixel_lor_hits_ptr_;
  }

  void compact_pixel_index(S i_pixel) {
    if (!pixel_lor_hits_ptr_[i_pixel] || i_pixel >= n_pixels_)
      return;

    // ensure we have enough space for the all LORs for that pixel
    pixel_lor_hits_[i_pixel].resize(pixel_lor_count_[i_pixel]);

    for (S i_lor = 0, lor_count = 0; i_lor < n_lors_; ++i_lor) {
      for (S position = 0; position < this->n_tof_positions(); ++position) {
        auto hits = pixel_lor_hits_ptr_[i_pixel][
            i_lor * this->n_tof_positions() + position];
        if (hits > 0) {
          LOR lor = index_to_lor_[i_lor];
          if (::abs(lor.first - lor.second) < 1) {
            std::ostringstream msg;
            msg << __PRETTY_FUNCTION__ << " invalid LOR " << i_lor << " ("
                << lor.first << ", " << lor.second << ") for pixel index "
                << i_pixel;
            throw(msg.str());
          }
          pixel_lor_hits_[i_pixel][lor_count++] = SparseElement(
              index_to_lor_[i_lor], position, index_to_pixel_[i_pixel], hits);
        }
      }
    }

    delete[] pixel_lor_hits_ptr_[i_pixel], pixel_lor_hits_ptr_[i_pixel] = NULL;

    std::sort(pixel_lor_hits_[i_pixel].begin(),
              pixel_lor_hits_[i_pixel].end(),
              SparseElementLORComparator());
  }

  SparseMatrix to_sparse() {
    SparseMatrix sparse(this->n_pixels_in_row(),
                        this->n_detectors(),
                        this->n_emissions(),
                        true,
                        this->n_tof_positions() > 0);
    sparse.reserve(size_);
    for (S i_pixel = 0; i_pixel < n_pixels_; ++i_pixel) {
      for (auto it = pixel_lor_hits_[i_pixel].begin();
           it != pixel_lor_hits_[i_pixel].end();
           ++it) {
        sparse.push_back(*it);
      }
    }
    return sparse;
  }

  MatrixPixelMajor& operator<<(SparseMatrix& sparse) {

    sparse.sort_by_pixel();
    this->add_emissions(sparse.n_emissions());

    S i_current_pixel = n_pixels_;

    for (auto it = sparse.begin(); it != sparse.end(); ++it) {
      Pixel pixel = it->pixel;
      if (!sparse.triangular()) {
        pixel.x -= n_pixels_in_row_half_;
        pixel.y -= n_pixels_in_row_half_;
        if (pixel.x < 0 || pixel.y < 0 || pixel.y < pixel.x)
          continue;
      }
      auto lor = it->lor;
      auto hits = it->hits;
      S i_pixel = pixel.index();

      if (i_current_pixel != i_pixel) {
        compact_pixel_index(i_current_pixel);
        i_current_pixel = i_pixel;
      }
      hit_lor(lor, it->position, i_pixel, hits);
      this->hit(i_pixel, hits);
    }
    compact_pixel_index(i_current_pixel);
    return *this;
  }

  // for testing purposes
  S lor_hits_at_pixel_index(LOR lor, S i_pixel) {
    auto it =
        std::lower_bound(pixel_lor_hits_[i_pixel].begin(),
                         pixel_lor_hits_[i_pixel].end(),
                         SparseElement(lor, 0, index_to_pixel_[i_pixel], 0),
                         SparseElementLORComparator());

    if (it == pixel_lor_hits_[i_pixel].end())
      return 0;
    return it->hits;
  }
  S size() const { return size_; }
  S n_lors_at_pixel_index(S i_pixel) const { return pixel_lor_count_[i_pixel]; }
  S n_pixels() const { return n_pixels_; }

 private:
  // disable copy contructor
  MatrixPixelMajor(const MatrixPixelMajor& rhs __attribute__((unused)))
      : Super(0, 0) {
    throw(__PRETTY_FUNCTION__);
  }

  struct SparseElementLORComparator {
    bool operator()(const SparseElement& a, const SparseElement& b) const {
      return a.lor < b.lor;
    }
  };

  S n_pixels_in_row_half_;
  S n_pixels_;
  S n_lors_;
  S size_;
  Hit** pixel_lor_hits_ptr_;
  std::vector<std::vector<SparseElement>> pixel_lor_hits_;
  std::vector<S> pixel_lor_count_;
  std::vector<LOR> index_to_lor_;
  std::vector<Pixel> index_to_pixel_;
};
