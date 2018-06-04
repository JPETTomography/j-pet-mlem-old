#pragma once

#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <type_traits>

#include "lor.h"
#include "2d/geometry/pixel.h"

#if !__CUDACC__
#include "ring_scanner.h"
#include "sparse_matrix.h"
#include "2d/geometry/pixel_map.h"
#endif

// This makes every iteration calculate rho-detected,
// aka rho*sensitivity
#define USE_RHO_DETECTED 0

namespace PET2D {
namespace Barrel {

/// 2D barrel PET reconstruction
template <typename FType, typename SType, typename HitType>
class Reconstruction {
 public:
  using F = FType;
  using S = SType;
  using Size = typename std::common_type<S, int>::type;
  using Hit = HitType;
  using Pixel = PET2D::Pixel<S>;
  using LOR = Barrel::LOR<S>;
  /// Input element for reconstruction
  typedef struct {
    LOR lor;
    S position;
    Hit mean;
  } Mean;

#if !__CUDACC__
  using Means = std::vector<Mean>;
  using Output = PET2D::PixelMap<Pixel, F>;
  using Sensitivity = PET2D::PixelMap<Pixel, Hit>;
  using Matrix = SparseMatrix<Pixel, LOR, Hit>;

  Reconstruction(Matrix& matrix, bool use_sensitivity = true)
      : n_pixels_in_row_(matrix.n_pixels_in_row()),
        total_n_pixels_(n_pixels_in_row_ * n_pixels_in_row_),
        sensitivity_(n_pixels_in_row_, n_pixels_in_row_),
        scale_(n_pixels_in_row_, n_pixels_in_row_, 1),
        rho_(n_pixels_in_row_, n_pixels_in_row_, 1),
        matrix_(matrix) {

    if (use_sensitivity) {
      for (const auto element : matrix_) {
        sensitivity_[element.pixel] += element.hits;
      }

      auto n_emissions = static_cast<F>(matrix.n_emissions());

      for (Size p = 0; p < total_n_pixels_; ++p) {
        Hit pixel_sensitivity = sensitivity_[p];
        if (pixel_sensitivity > 0) {
          scale_[p] = static_cast<F>(n_emissions) / pixel_sensitivity;
        }
      }
    }

    matrix_.sort_by_lor();
  }

  /// Reads means from given input stream
  Reconstruction& operator<<(std::istream& in_means) {
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

    if (matrix_.n_tof_positions() > 1) {
      std::sort(means_.begin(), means_.end(), SortByLORNPosition());
    } else {
      std::sort(means_.begin(), means_.end(), SortByLOR());
    }

    return *this;
  }

  /// Performs n_iterations of the EMT algorithm
  template <typename ProgressCallback>
  void operator()(ProgressCallback& progress,  ///< progress callback
                  int n_iterations,            ///< iterations to perform
                  int n_iterations_so_far = 0  ///< iterations so far
                  ) {
    F* y = (F*)alloca(total_n_pixels_ * sizeof(F));

#if USE_RHO_DETECTED
    // tranform rho -> rho-detected
    // starting from constant 1 makes effectively rho=sensitivity
    for (Size p = 0; p < total_n_pixels_; ++p) {
      rho_[p] /= scale_[p];
    }
#endif

    for (Size iteration = 0; iteration < n_iterations; ++iteration) {
      progress(iteration + n_iterations_so_far);

      for (Size p = 0; p < total_n_pixels_; ++p) {
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

        // skip LORs & positions that does not exist in system matrix
        while (means_it != means_.end() &&
               (matrix_it->lor > means_it->lor ||
                (matrix_it->lor == means_it->lor &&
                 matrix_it->position > means_it->position))) {
          if (iteration == 0) {
            // this warning should not appear if system matrix is complete
            std::cerr << "warning: mean LOR (" << means_it->lor.first << ", "
                      << means_it->lor.second << ") position "
                      << means_it->position << " not found in system matrix"
                      << std::endl;
          }
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
            auto p = pixel_index(matrix_it->pixel);
            u += rho_[p] * static_cast<F>(matrix_it->hits)
#if USE_RHO_DETECTED
                 * scale_[p];
#endif
            ;
            ++matrix_it;
          }
          F phi = means_it->mean / u;

          // count y for current lor
          matrix_it = prev_it;
          while (matrix_it->lor == lor && matrix_it->position == position) {
            auto p = pixel_index(matrix_it->pixel);
            y[p] += phi * static_cast<F>(matrix_it->hits)
#if USE_RHO_DETECTED
                    * scale_[p];
#endif
            ;
            ++matrix_it;
          }
        } else {
          // skip this LOR
          while (matrix_it->lor == lor && matrix_it->position == position)
            ++matrix_it;
        }
        ++means_it;
      }

      for (Size p = 0; p < total_n_pixels_; ++p) {
        rho_[p] *= y[p]
#if !USE_RHO_DETECTED
                   * scale_[p]
#endif
            ;
      }

      progress(iteration + n_iterations_so_far, true);
    }

#if USE_RHO_DETECTED
    // apply sensitivity for final rho-detected -> rho
    for (Size p = 0; p < total_n_pixels_; ++p) {
      rho_[p] *= scale_[p];
    }
#endif
  }

  S n_pixels_in_row() { return n_pixels_in_row_; }
  F rho(const S p) const { return rho_[p]; }
  F rho(const Pixel& pixel) const { return rho_[pixel_index(pixel)]; }
  const Sensitivity& sensitivity() const { return sensitivity_; }
  const Output& rho() const { return rho_; }
  Output rho_detected() const {
    Output rho_detected(rho_);
    for (Size p = 0; p < total_n_pixels_; ++p) {
      rho_detected[p] /= scale_[p];
    }
    return rho_detected;
  }
  const Output& scale() const { return scale_; }
  const Means& means() const { return means_; }

 private:
  Size pixel_index(const Pixel& p) const {
    return p.y * static_cast<Size>(n_pixels_in_row_) + p.x;
  }

  S n_pixels_in_row_;        ///< number of pixels in row in image
  Size total_n_pixels_;      ///< total number of pixel in image
  Sensitivity sensitivity_;  ///< pixel sensitivity
  Output scale_;             ///< inverse of sensitivity
  Output rho_;               ///< reconstruction output
  Matrix& matrix_;           ///< system matrix used for reconstruction
  Means means_;              ///< input means (eg. from phantom simmulation)

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
#endif
};
}  // Barrel
}  // PET2D
