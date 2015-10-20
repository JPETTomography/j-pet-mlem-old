#pragma once

#include "lor.h"
#include "../geometry/pixel.h"
#if !__CUDACC__
#include "lor_geometry.h"
#include "sparse_matrix.h"
#include "geometry.h"
#endif

namespace PET2D {
namespace Barrel {

/// Keeps flat information about barrel, including grid and all LOR info
////
/// This class is simplification of PET2D::Barrel::Geometry, the difference
/// is that it keeps geometry as flat array of \c PixelInfo structures, not
/// using \c std::vector or any STL containers.
///
/// \see PET2D::Barrel::Geometry
/// \see PET2D::Barrel::SparseMatrix
template <typename FType, typename SType, typename HitType>
class SimpleGeometry {
 public:
  using F = FType;
  using S = SType;
  using Hit = HitType;
  using Pixel = PET2D::Pixel<S>;
  using LOR = PET2D::Barrel::LOR<S>;
#if !__CUDACC__
  using SparseMatrix = PET2D::Barrel::SparseMatrix<Pixel, LOR, Hit>;
  using Geometry = PET2D::Barrel::Geometry<F, S>;
#endif

  /// Information about single pixel in context of certain LOR
  struct PixelInfo {
    Pixel pixel;  ///< pixel itself
    F t;          ///< position of projected Pixel along the LOR
    F weight;     ///< custom weight of the Pixel eg. probability
  };

  /// Construct empty geometry information for given number of detectors.
  SimpleGeometry(S n_detectors,        ///< number of detectors
                 size_t n_pixel_infos  ///< total number of pixel infos
                 )
      : n_detectors(n_detectors),
        n_lors(((size_t(n_detectors) + 1) * size_t(n_detectors)) / 2),
        n_pixel_infos(n_pixel_infos),
        pixel_infos(new PixelInfo[n_pixel_infos]),
        lor_pixel_info_start(new size_t[n_lors]),
        lor_pixel_info_end(new size_t[n_lors]) {}

  ~SimpleGeometry() {
    delete[] pixel_infos;
    delete[] lor_pixel_info_start;
    delete[] lor_pixel_info_end;
  }

#if !__CUDACC__
  /// Construct geometry information out of sparse matrix.
  SimpleGeometry(const SparseMatrix& sparse_matrix  ///< sparse matrix
                 )
      : SimpleGeometry(sparse_matrix.n_detectors(), sparse_matrix.size()) {
    const auto end_lor = LOR::end_for_detectors(n_detectors);
    auto lor = end_lor;
    auto lor_index = lor.index();
    if (sparse_matrix.sorted() != SparseMatrix::BY_LOR_N_PIXEL)
      throw("input sparse matrix must be sorted by lor and pixel");
    size_t index = 0;
    for (const auto& element : sparse_matrix) {
      // check if we have new LOR
      if (element.lor != lor) {
        if (lor != end_lor) {
          lor_pixel_info_end[lor_index] = index;
        }
        lor = element.lor;
        lor_index = lor.index();
        lor_pixel_info_start[lor_index] = index;
      }
      // assign information for this pixel info
      pixel_infos[index].pixel = element.pixel;
      pixel_infos[index].weight =
          (F)element.hits / (F)sparse_matrix.n_emissions();
      ++index;
    }
  }

  /// Construct geometry information out of more advanced geometry.
  ////
  /// Takes PET2D::Barrel::Geometry class instance and flattens it.
  SimpleGeometry(const Geometry& geometry  ///< advanced geometry
                 )
      : SimpleGeometry(geometry.n_detectors, geometry.n_pixel_infos()) {
    size_t size = 0;
    for (const auto& lor_geometry : geometry) {
      const auto& lor = lor_geometry.lor;
      auto lor_index = lor.index();
      lor_pixel_info_start[lor_index] = size;
      for (const auto& geometry_pixel_info : lor_geometry.pixel_infos) {
        auto& pixel_info = pixel_infos[size++];
        pixel_info.pixel = geometry_pixel_info.pixel;
        pixel_info.t = geometry_pixel_info.t;
        pixel_info.weight = geometry_pixel_info.weight;
      }
      lor_pixel_info_end[lor_index] = size;
    }
  }
#endif

  const S n_detectors;                 ///< number of detectors
  const size_t n_lors;                 ///< total number of LORs
  const size_t n_pixel_infos;          ///< total number of pixel infos
  PixelInfo* const pixel_infos;        ///< pointer to pixel infos array
  size_t* const lor_pixel_info_start;  ///< pointer to array holding start of
                                       ///  pixel infos for given LOR
  size_t* const lor_pixel_info_end;    ///< pointer to array holding end of
                                       ///  pixel infos for given LOR
};

}  // Barrel
}  // PET2D
