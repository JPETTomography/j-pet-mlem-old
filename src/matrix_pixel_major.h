#pragma once

#include "triangular_pixel_map.h"

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

template <typename LORType, typename FType = double, typename SType = int>
class MatrixPixelMajor : public TriangularPixelMap<FType, SType> {
 public:
  typedef LORType LOR;
  typedef FType F;
  typedef SType S;
  typedef typename std::make_signed<S>::type SS;
  typedef TriangularPixelMap<F, S> Super;
  typedef std::pair<LOR, S> Hit;
  typedef std::pair<std::pair<LOR, S>, S> Pair;

  MatrixPixelMajor(S n_pixels, S n_detectors)
      : TriangularPixelMap<F, S>(n_pixels),
        end_(LOR::end_for_detectors(n_detectors)),
        n_pixels_(n_pixels),
        n_pixels_half_(n_pixels_ / 2),
        total_n_pixels_(n_pixels_half_ * (n_pixels_half_ + 1) / 2),
        n_lors_(LOR::end_for_detectors(n_detectors).t_index()),
        n_entries_(0),
        pixel_tmp_(total_n_pixels_, NULL),
        pixel_(total_n_pixels_),
        pixel_count_(total_n_pixels_),
        index_to_lor_(LOR::end_for_detectors(n_detectors).t_index()) {
    for (auto lor = begin(); lor != end(); ++lor) {
      index_to_lor_[lor.t_index()] = lor;
    }
  }

  static LOR begin() { return LOR(); }
  const LOR end() { return end_; }

  S n_entries() const { return n_entries_; }
  S n_lors(S p) const { return pixel_count_[p]; }
  S total_n_pixels() const { return total_n_pixels_; }

  void add_to_t_matrix(const LOR& lor, S i_pixel) {

    if (!pixel_tmp_[i_pixel]) {
      pixel_tmp_[i_pixel] = new S[n_lors_]();
    }
    if (pixel_tmp_[i_pixel][t_lor_index(lor)] == 0) {
      pixel_count_[i_pixel]++;
      n_entries_++;
    }
    pixel_tmp_[i_pixel][t_lor_index(lor)]++;
  }

  ~MatrixPixelMajor() {
    for (S p = 0; p < total_n_pixels_; ++p) {
      if (pixel_tmp_[p]) {
        delete[] pixel_tmp_[p];
      }
    }
  }

  S t_get_element(LOR lor, S i_pixel) {
    auto hit = std::lower_bound(pixel_[i_pixel].begin(),
                                pixel_[i_pixel].end(),
                                std::make_pair(lor, 0),
                                HitComparator());

    if (hit == pixel_[i_pixel].end())
      return 0;
    return (*hit).second;
  }

  void finalize_pixel(S i_pixel) {
    pixel_[i_pixel].resize(pixel_count_[i_pixel]);
    S it = 0;
    for (S lor = 0; lor < n_lors_; ++lor) {
      S hits;
      if ((hits = pixel_tmp_[i_pixel][lor]) > 0) {
        pixel_[i_pixel][it] = std::make_pair(index_to_lor_[lor], hits);
        it++;
      }
    }
    delete[] pixel_tmp_[i_pixel];
    pixel_tmp_[i_pixel] = NULL;
    std::sort(pixel_[i_pixel].begin(), pixel_[i_pixel].end(), HitComparator());
  }

  static S t_lor_index(const LOR& lor) { return lor.t_index(); }

  void to_pairs() {
    pair_.reserve(n_entries());
    for (S p = 0; p < total_n_pixels(); ++p) {
      for (auto it = pixel_[p].begin(); it != pixel_[p].end(); ++it) {
        pair_.push_back(
            std::make_pair(std::make_pair((*it).first, p), (*it).second));
      }
    }
  }

  Pair pair(S p) const { return pair_[p]; }

  void sort_pairs_by_lors() {
    std::sort(pair_.begin(), pair_.end(), LorSorter());
  }

  void sort_pairs_by_pixels() {
    std::sort(pair_.begin(), pair_.end(), PixelSorter());
  }

 private:
  struct HitComparator {
    bool operator()(const Hit& a, const Hit& b) const {
      return a.first < b.first;
    }
  };

  struct LorSorter {
    bool operator()(const Pair& a, const Pair& b) const {
      return a.first.first < b.first.first;
    }
  };

  struct PixelSorter {
    bool operator()(const Pair& a, const Pair& b) const {
      return a.first.second < b.first.second;
    }
  };

  LOR end_;
  S n_pixels_;
  S n_pixels_half_;
  S total_n_pixels_;
  S n_lors_;
  S n_entries_;
  std::vector<S*> pixel_tmp_;
  std::vector<std::vector<Hit>> pixel_;
  std::vector<S> pixel_count_;
  std::vector<LOR> index_to_lor_;
  std::vector<Pair> pair_;
};
