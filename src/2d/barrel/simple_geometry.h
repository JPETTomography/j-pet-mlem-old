#pragma once

#include "lor.h"
#include "../geometry/pixel.h"
#include "../geometry/line_segment.h"
#if !__CUDACC__
#include "lor_geometry.h"
#include "sparse_matrix.h"
#include "geometry.h"
#endif

namespace PET2D {
namespace Barrel {

/// Keeps flat SOA information about barrel, including grid and all LOR info
////
/// This class is simplification of PET2D::Barrel::Geometry, the difference
/// is that it keeps geometry as flat structure of arrays equivalent of
/// \c PixelInfo structures, not using \c std::vector or any STL containers.
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
  using LineSegment = PET2D::LineSegment<F>;
#if !__CUDACC__
  using SparseMatrix = PET2D::Barrel::SparseMatrix<Pixel, LOR, Hit>;
  using Geometry = PET2D::Barrel::Geometry<F, S>;
#endif

  /// Construct empty geometry information for given number of detectors.
  SimpleGeometry(S n_detectors,        ///< number of detectors
                 size_t n_pixel_infos  ///< total number of pixel infos
                 )
      : n_detectors(n_detectors),
        n_lors(((size_t(n_detectors) + 1) * size_t(n_detectors)) / 2),
        lor_line_segments(new LineSegment[n_lors]),
        n_pixel_infos(n_pixel_infos),
        pixels(new Pixel[n_pixel_infos]),
        pixel_positions(new F[n_pixel_infos]),
        pixel_weights(new F[n_pixel_infos]),
        lor_pixel_info_begin(new size_t[n_lors]),
        lor_pixel_info_end(new size_t[n_lors]) {}

  ~SimpleGeometry() {
    delete[] lor_line_segments;
    delete[] pixels;
    delete[] pixel_positions;
    delete[] pixel_weights;
    delete[] lor_pixel_info_begin;
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
        lor_pixel_info_begin[lor_index] = index;
      }
      // assign information for this pixel info
      pixels[index] = element.pixel;
      pixel_weights[index] = (F)element.hits / (F)sparse_matrix.n_emissions();
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
      lor_line_segments[lor_index] = lor_geometry.segment;
      lor_pixel_info_begin[lor_index] = size;
      for (const auto& geometry_pixel_info : lor_geometry.pixel_infos) {
        pixels[size] = geometry_pixel_info.pixel;
        pixel_positions[size] = geometry_pixel_info.t;
        pixel_weights[size] = geometry_pixel_info.weight;
        ++size;
      }
      lor_pixel_info_end[lor_index] = size;
    }
  }
#endif

  const S n_detectors;                   ///< number of detectors
  const size_t n_lors;                   ///< total number of LORs
  LineSegment* const lor_line_segments;  ///< LOR line segments array
  const size_t n_pixel_infos;            ///< total number of pixel infos
  Pixel* const pixels;                   ///< pointer to pixels array
  F* const pixel_positions;              ///< projected pixels along the LOR
  F* const pixel_weights;                ///< pixel weights eg. probability
  size_t* const lor_pixel_info_begin;    ///< pointer to array holding start of
                                         ///  pixel infos for given LOR
  size_t* const lor_pixel_info_end;      ///< pointer to array holding end of
                                         ///  pixel infos for given LOR
};

}  // Barrel
}  // PET2D
