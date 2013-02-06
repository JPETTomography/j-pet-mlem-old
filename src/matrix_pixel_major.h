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
/// It alos means that as different threads are processing different pixels,
/// there is no need for synchronisation in add_pixels.
///
/// The reconstruction is however done using the lor_major matrix.
/// So this class has the possibility to write the matrix down in the
/// triangular lor_major and full lor_major format.

template <typename LORType, typename SType = int, typename HitType = int>
class MatrixPixelMajor : public Matrix<LORType, SType, HitType> {
  typedef Matrix<LORType, SType, HitType> Super;

 public:
  typedef LORType LOR;
  typedef SType S;
  typedef HitType Hit;
  typedef typename std::make_signed<S>::type SS;
  typedef std::pair<LOR, Hit> LORHit;
  typedef std::tuple<LOR, S, Hit> LORPixelHit;
  typedef std::vector<LOR> LORList;
  typedef std::vector<LORPixelHit> LORPixelHitList;
  typedef typename Super::Sparse Sparse;

  MatrixPixelMajor(S n_pixels, S n_detectors)
      : Super(n_pixels, n_detectors),
        n_pixels_(n_pixels),
        n_pixels_half_(n_pixels_ / 2),
        total_n_pixels_(n_pixels_half_ * (n_pixels_half_ + 1) / 2),
        n_lors_(LOR::end_for_detectors(n_detectors).t_index()),
        size_(0),
        pixel_lor_hits_ptr_(new S* [total_n_pixels_]()),
        pixel_lor_hits_(total_n_pixels_),
        pixel_lor_count_(total_n_pixels_),
        index_to_lor_(LOR::end_for_detectors(n_detectors).t_index()) {

    // store index to LOR mapping
    for (auto lor = Super::begin(); lor != Super::end(); ++lor) {
      index_to_lor_[lor.t_index()] = lor;
    }
  }

  S size() const { return size_; }

  S n_lors(S p) const { return pixel_lor_count_[p]; }
  S total_n_pixels() const { return total_n_pixels_; }

  void add_to_t_matrix(const LOR& lor, S i_pixel) {

    if (!pixel_lor_hits_ptr_[i_pixel]) {
      pixel_lor_hits_ptr_[i_pixel] = new S[n_lors_]();
    }
    if (pixel_lor_hits_ptr_[i_pixel][t_lor_index(lor)] == 0) {
      pixel_lor_count_[i_pixel]++;
      size_++;
    }
    pixel_lor_hits_ptr_[i_pixel][t_lor_index(lor)]++;
  }

  ~MatrixPixelMajor() {
    for (S p = 0; p < total_n_pixels_; ++p) {
      if (pixel_lor_hits_ptr_[p]) {
        delete[] pixel_lor_hits_ptr_[p];
      }
    }
  }

  S t_get_element(LOR lor, S i_pixel) {
    auto it = std::lower_bound(pixel_lor_hits_[i_pixel].begin(),
                               pixel_lor_hits_[i_pixel].end(),
                               LORHit(lor, 0),
                               LORHitComparator());

    if (it == pixel_lor_hits_[i_pixel].end())
      return 0;
    return it->second;
  }

  void compact_pixel_index(S i_pixel) {
    // ensure we have enough space for the all LORs for that pixel
    pixel_lor_hits_[i_pixel].resize(pixel_lor_count_[i_pixel]);

    for (S i_lor = 0, lor_count = 0; i_lor < n_lors_; ++i_lor) {
      S hits;
      if ((hits = pixel_lor_hits_ptr_[i_pixel][i_lor]) > 0) {
        pixel_lor_hits_[i_pixel][lor_count++] =
            LORHit(index_to_lor_[i_lor], hits);
      }
    }

    delete[] pixel_lor_hits_ptr_[i_pixel], pixel_lor_hits_ptr_[i_pixel] = NULL;

    std::sort(pixel_lor_hits_[i_pixel].begin(),
              pixel_lor_hits_[i_pixel].end(),
              LORHitComparator());
  }

  static S t_lor_index(const LOR& lor) { return lor.t_index(); }

  Sparse to_sparse() {
    Sparse sparse;
    sparse.reserve(size_);
    for (S p = 0; p < total_n_pixels_; ++p) {
      for (auto it = pixel_lor_hits_[p].begin(); it != pixel_lor_hits_[p].end();
           ++it) {
        sparse.push_back(LORPixelHit(it->first, p, it->second));
      }
    }
    return sparse;
  }

 private:
  struct LORHitComparator {
    bool operator()(const LORHit& a, const LORHit& b) const {
      return a.first < b.first;
    }
  };

  S n_pixels_;
  S n_pixels_half_;
  S total_n_pixels_;
  S n_lors_;
  S size_;
  Hit** pixel_lor_hits_ptr_;
  std::vector<std::vector<LORHit>> pixel_lor_hits_;
  std::vector<S> pixel_lor_count_;
  LORList index_to_lor_;
};
