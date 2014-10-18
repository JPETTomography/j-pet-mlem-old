#pragma once

#include <iostream>
#include <vector>
#include <list>
#include <algorithm>

#include "detector_ring.h"
#include "sparse_matrix.h"
#include "2d/geometry/pixel.h"
#include "lor.h"

template <typename FType = double, typename SType = int> class Reconstruction {
 public:
  typedef FType F;
  typedef SType S;
  typedef ::Pixel<S> Pixel;
  typedef ::LOR<S> LOR;
  typedef struct {
    LOR lor;
    S position;
    S mean;
  } Mean;
  typedef std::vector<Mean> Means;
  typedef std::vector<F> Output;
  typedef SparseMatrix<Pixel, LOR> Matrix;

  Reconstruction(S n_iterations,
                 Matrix& matrix,
                 std::istream& in_means,
                 F threshold = static_cast<F>(0))
      : threshold(threshold), n_iterations_(n_iterations), matrix_(matrix) {

    n_pixels_in_row_ = matrix_.n_pixels_in_row();
    n_emissions_ = matrix_.n_emissions();
    n_detectors_ = matrix_.n_detectors();

    total_n_pixels_ = n_pixels_in_row_ * n_pixels_in_row_;

    rho_.resize(total_n_pixels_);
    rho_detected_.resize(total_n_pixels_, static_cast<F>(1));
    scale_.resize(total_n_pixels_, static_cast<F>(0));

    for (auto it = matrix_.begin(); it != matrix_.end(); ++it) {
      scale_[pixel_index(it->pixel)] += it->hits;
    }

    auto n_emissions = static_cast<F>(n_emissions_);

    for (S p = 0; p < total_n_pixels_; ++p) {
      if (scale_[p] > 0) {
        scale_[p] = n_emissions / scale_[p];
      }
    }

    for (;;) {
      Mean mean;
      in_means >> mean.lor.first >> mean.lor.second;
      if (mean.lor.first < mean.lor.second)
        std::swap(mean.lor.first, mean.lor.second);
      if (matrix_.n_tof_positions() > 1) {
        in_means >> mean.position;
      } else {
        mean.position = 0;
      }
      in_means >> mean.mean;
      if (in_means.eof())
        break;
      means_.push_back(mean);
    }

    matrix_.sort_by_lor();

    if (matrix_.n_tof_positions() > 1) {
      std::sort(means_.begin(), means_.end(), SortByLORNPosition());
    } else {
      std::sort(means_.begin(), means_.end(), SortByLOR());
    }
  }

  void emt(S n_iterations) {
    F* y = (F*)alloca(n_pixels_in_row_ * n_pixels_in_row_);

    for (S i = 0; i < n_iterations; ++i) {
      std::cout << ".", std::cout.flush();

      for (S p = 0; p < total_n_pixels_; ++p) {
        y[p] = static_cast<F>(0);
      }

      auto matrix_it = matrix_.begin();
      auto means_it = means_.begin();
      for (;;) {
        // skip LORs that does not exist in means
        while (matrix_it != matrix_.end() &&
               (matrix_it->lor < means_it->lor ||
                (matrix_it->lor == means_it->lor &&
                 matrix_it->position < means_it->position))) {
#if DEBUG
          std::cerr << "skip matrix LOR (" << matrix_it->lor.first << ", "
                    << matrix_it->lor.second << ") position "
                    << matrix_it->position << std::endl;
#endif
          ++matrix_it;
        }

        // skip LORs that does not exist in system matrix
        while (means_it != means_.end() &&
               (matrix_it->lor > means_it->lor ||
                (matrix_it->lor == means_it->lor &&
                 matrix_it->position > means_it->position))) {
#if 1  // this warning should not appear if system matrix is complete
          std::cerr << "warning: mean LOR (" << means_it->lor.first << ", "
                    << means_it->lor.second << ") position "
                    << means_it->position << " not found in system matrix"
                    << std::endl;
#endif
          ++means_it;
        }

        // check if we are EOT
        if (matrix_it == matrix_.end() || means_it == means_.end())
          break;

        if (matrix_it->lor != means_it->lor ||
            matrix_it->position != means_it->position)
          continue;

        // store current lor & position
        auto lor = matrix_it->lor;
        auto position = matrix_it->position;

#if DEBUG
        std::cerr << "processing LOR (" << lor.first << ", " << lor.second
                  << ") position " << position << std::endl;
#endif

        // if there any mean anyway here?
        if (means_it->mean > 0) {
          F u = static_cast<F>(0);
          auto prev_it = matrix_it;

          // count u for current LOR
          while (matrix_it->lor == lor && matrix_it->position == position) {
            auto i_pixel = pixel_index(matrix_it->pixel);
            u += rho_detected_[i_pixel] * static_cast<F>(matrix_it->hits) *
                 scale_[i_pixel];
            ++matrix_it;
          }
          F phi = means_it->mean / u;

          // count y for current lor
          matrix_it = prev_it;
          while (matrix_it->lor == lor && matrix_it->position == position) {
            auto i_pixel = pixel_index(matrix_it->pixel);
            y[pixel_index(matrix_it->pixel)] +=
                phi * static_cast<F>(matrix_it->hits) * scale_[i_pixel];
            ++matrix_it;
          }
        } else {
          // skip this LOR
          while (matrix_it->lor == lor && matrix_it->position == position)
            ++matrix_it;
        }
        ++means_it;
      }

      for (S p = 0; p < total_n_pixels_; ++p) {
        rho_detected_[p] *= y[p];
      }
    }

    for (S p = 0; p < total_n_pixels_; ++p) {
      rho_[p] = rho_detected_[p] * scale_[p];
    }
  }

  S n_pixels_in_row() { return n_pixels_in_row_; }
  F rho(const S p) const { return rho_[p]; }
  F rho(const Pixel& pixel) const { return rho_[pixel_index(pixel)]; }
  Output rho() const { return rho_; }
  Output rho_detected() { return rho_detected_; }

 public:
  F threshold;

 private:
  S pixel_index(const Pixel& p) const { return p.y * n_pixels_in_row_ + p.x; }

  S n_detectors_;
  S n_pixels_in_row_;
  S total_n_pixels_;
  S n_iterations_;
  S n_emissions_;
  Output scale_;
  Output rho_;
  Output rho_detected_;
  Matrix& matrix_;
  Means means_;

  struct SortByLOR {
    bool operator()(const Mean& a, const Mean& b) const {
      return a.lor < b.lor;
    }
  };

  struct SortByLORNPosition {
    bool operator()(const Mean& a, const Mean& b) const {
      return a.lor < b.lor || (a.lor == b.lor && a.position < b.position);
    }
  };
};
