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

  struct PixelInfo {
    Pixel pixel;  ///< pixel itself
    F t;          ///< position of projected Pixel along the LOR
    F weight;     ///< custom weight of the Pixel eg. probability
  };

  SimpleGeometry(S n_detectors, size_t n_pixel_infos)
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
  SimpleGeometry(const SparseMatrix& sparse_matrix)
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

  SimpleGeometry(const Geometry& geometry)
      : SimpleGeometry(geometry.n_detectors(), geometry.n_pixel_infos()) {
    size_t size = 0;
    for (const auto& lor_geometry : geometry) {
      const auto& lor = lor_geometry.lor;
      auto lor_index = lor.index();
      lor_pixel_info_start[lor_index] = size;
      for (const auto& geometry_pixel_info : lor_geometry) {
        auto& pixel_info = pixel_infos[size++];
        pixel_info.pixel = geometry_pixel_info.pixel;
        pixel_info.t = geometry_pixel_info.t;
        pixel_info.weight = geometry_pixel_info.weight;
      }
      lor_pixel_info_end[lor_index] = size;
    }
  }
#endif

  const S n_detectors;
  const size_t n_lors;
  const size_t n_pixel_infos;
  PixelInfo* const pixel_infos;
  size_t* const lor_pixel_info_start;
  size_t* const lor_pixel_info_end;
};

}  // Barrel
}  // PET2D
