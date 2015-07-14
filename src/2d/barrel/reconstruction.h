#pragma once

#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <type_traits>

#include "ring_scanner.h"
#include "sparse_matrix.h"
#include "2d/geometry/pixel.h"
#include "lor.h"

#define DEBUG 0

namespace PET2D {
namespace Barrel {

/// 2D barrel PET reconstruction
template <typename FType, typename SType, typename HitType>
class Reconstruction {
 public:
  using F = FType;
  using S = SType;
  using I = typename std::common_type<S, int>::type;
  using Hit = HitType;
  using Pixel = PET2D::Pixel<S>;
  using LOR = Barrel::LOR<S>;
  /// Input element for reconstruction
  typedef struct {
    LOR lor;
    S position;
    Hit mean;
  } Mean;
  using Means = std::vector<Mean>;
  using Output = std::vector<F>;
  using Matrix = SparseMatrix<Pixel, LOR, S, Hit>;

  Reconstruction(Matrix& matrix,
                 std::istream& in_means,
                 bool use_sensitivity = true)
      : matrix_(matrix) {

    n_pixels_in_row_ = matrix_.n_pixels_in_row();
    n_emissions_ = matrix_.n_emissions();
    n_detectors_ = matrix_.n_detectors();

    total_n_pixels_ = n_pixels_in_row_ * n_pixels_in_row_;

    rho_.resize(total_n_pixels_);
    rho_detected_.resize(total_n_pixels_, 1);

    if (use_sensitivity) {
      std::vector<Hit> sensitivity(total_n_pixels_, 0);
      scale_.resize(total_n_pixels_, 0);

      for (const auto element : matrix_) {
        sensitivity[pixel_index(element.pixel)] += element.hits;
      }

      auto n_emissions = static_cast<F>(n_emissions_);

      for (I p = 0; p < total_n_pixels_; ++p) {
        Hit pixel_sensitivity = sensitivity[p];
        if (pixel_sensitivity > 0) {
          scale_[p] = static_cast<F>(n_emissions) / pixel_sensitivity;
        }
      }
    } else {
      scale_.resize(total_n_pixels_, 1);
    }

    // Read the mean (detector response file)
    for (;;) {
      Mean mean;
      in_means >> mean.lor.first >> mean.lor.second >> mean.position >>
          mean.mean;
      if (in_means.eof())
        break;
      if (mean.lor.first < mean.lor.second)
        std::swap(mean.lor.first, mean.lor.second);
      means_.push_back(mean);
    }

    matrix_.sort_by_lor();

    if (matrix_.n_tof_positions() > 1) {
      std::sort(means_.begin(), means_.end(), SortByLORNPosition());
    } else {
      std::sort(means_.begin(), means_.end(), SortByLOR());
    }
  }

  void emt(I n_iterations) {
    F* y = (F*)alloca(n_pixels_in_row_ * n_pixels_in_row_ * sizeof(F));

    for (I i = 0; i < n_iterations; ++i) {
      std::cout << ".", std::cout.flush();

      for (I p = 0; p < total_n_pixels_; ++p) {
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
          ++matrix_it;
        }

        // skip LORs that does not exist in system matrix
        while (means_it != means_.end() && (matrix_it->lor > means_it->lor)) {
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

        // if there any mean anyway here?
        if (means_it->mean > 0) {
          F u = static_cast<F>(0);
          auto prev_it = matrix_it;

          // count u for current LOR
          while (matrix_it->lor == lor && matrix_it->position == position) {
            auto i_pixel = pixel_index(matrix_it->pixel);
            u += rho_detected_[i_pixel] * static_cast<F>(matrix_it->hits)  //
                 * scale_[i_pixel];
            ++matrix_it;
          }
          F phi = means_it->mean / u;

          // count y for current lor
          matrix_it = prev_it;
          while (matrix_it->lor == lor && matrix_it->position == position) {
            auto i_pixel = pixel_index(matrix_it->pixel);
            y[pixel_index(matrix_it->pixel)] +=
                phi * static_cast<F>(matrix_it->hits)  //
                * scale_[i_pixel];
            ++matrix_it;
          }
        } else {
          // skip this LOR
          while (matrix_it->lor == lor && matrix_it->position == position)
            ++matrix_it;
        }
        ++means_it;
      }

      for (I p = 0; p < total_n_pixels_; ++p) {
        rho_detected_[p] *= y[p];
      }
    }

    for (I p = 0; p < total_n_pixels_; ++p) {
      rho_[p] = rho_detected_[p] * scale_[p];
    }
  }

  S n_pixels_in_row() { return n_pixels_in_row_; }
  F rho(const S p) const { return rho_[p]; }
  F rho(const Pixel& pixel) const { return rho_[pixel_index(pixel)]; }
  const Output& rho() const { return rho_; }
  const Output& rho_detected() const { return rho_detected_; }
  const Output& scale() const { return scale_; }

 private:
  I pixel_index(const Pixel& p) const {
    return p.y * static_cast<I>(n_pixels_in_row_) + p.x;
  }

  S n_detectors_;
  S n_pixels_in_row_;
  I total_n_pixels_;
  I n_emissions_;
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
}  // Barrel
}  // PET2D
