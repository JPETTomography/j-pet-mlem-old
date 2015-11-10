#pragma once

#if !__CUDACC__
#include "util/bstream.h"
#include "util/read.h"
#endif

#include "lor_geometry.h"

namespace PET2D {
namespace Barrel {

/// \cond PRIVATE
#if !__CUDACC__
namespace {
template <typename FType> constexpr uint32_t magic() { return 0; }
template <> constexpr uint32_t magic<float>() { return "PETg"_4cc; }
template <> constexpr uint32_t magic<double>() { return "PETG"_4cc; }
}
#endif
/// \endcond

/// Keeps extended information about barrel, including grid and all LOR info
////
/// This class is complementary to PET2D::Barrel::SparseMatrix, the difference
/// is that it keeps geometry information for LOR-pixel pairs, like pixel
/// distance to LOR line, its position along LOR line, etc.
///
/// \see PET2D::Barrel::SparseMatrix
/// \see PET2D::Barrel::SimpleGeometry
template <typename FType, typename SType>
class Geometry : public std::vector<LORGeometry<FType, SType>> {
 public:
  using F = FType;
  using S = SType;

  using LORGeometry = PET2D::Barrel::LORGeometry<F, S>;
  using Base = std::vector<LORGeometry>;
  using Pixel = PET2D::Pixel<S>;
  using LOR = PET2D::Barrel::LOR<S>;
  using Grid = PET2D::PixelGrid<F, S>;
  using Point = PET2D::Point<F>;
  using PixelInfo = typename LORGeometry::PixelInfo;

  /// Construct geometry for given number of detector and grid description.
  Geometry(S n_detectors,    ///< total number of detectors
           const Grid& grid  ///< pixel grid description
           )
      : Base(((int(n_detectors) + 1) * (n_detectors)) / 2),
        n_detectors(n_detectors),
        grid(grid) {}

#if !__CUDACC__
  /// Construct geometry from stream
  Geometry(std::istream& in) : Geometry(util::read<uint32_t>(in), in) {}

  /// Construct geometry from binary stream
  Geometry(util::ibstream& in) : Geometry(in.read<uint32_t>(), in) {}

 private:
  // This private proxy constructor guarantee proper read order
  Geometry(S n_detectors, std::istream& in)
      : Base(((int(n_detectors) + 1) * (n_detectors)) / 2),
        n_detectors(n_detectors),
        grid(in) {
    for (;;) {
      LORGeometry lor_geometry(in);
      if (!in)
        break;
      (*this)[lor_geometry.lor] = std::move(lor_geometry);
    }
  }
  Geometry(uint32_t file_magic, S n_detectors, util::ibstream& in)
      : Base(((int(n_detectors) + 1) * (n_detectors)) / 2),
        n_detectors(n_detectors),
        grid(in) {
    if (file_magic != magic<F>())
      throw("unknown binary geometry file");
    for (;;) {
      LORGeometry lor_geometry(in);
      if (!in)
        break;
      (*this)[lor_geometry.lor] = std::move(lor_geometry);
    }
  }
  Geometry(uint32_t file_magic, util::ibstream& in)
      : Geometry(file_magic, in.read<S>(), in) {}

 public:
  /// Serialize geometry into output stream
  friend std::ostream& operator<<(std::ostream& out, const Geometry& geometry) {
    out << geometry.n_detectors << '\n' << geometry.grid;
    for (const auto& lor_geometry : geometry) {
      out << '\n' << lor_geometry;
    }
    return out;
  }

  /// Serialize geometry into binary output stream
  friend util::obstream& operator<<(util::obstream& out,
                                    const Geometry& geometry) {
    out << magic<F>();
    out << geometry.n_detectors;
    out << geometry.grid;
    for (const auto& lor_geometry : geometry) {
      out << lor_geometry;
    }
    return out;
  }
#endif

  /// Returns LOR geometry description for given LOR.
  LORGeometry& operator[](const LOR& lor) { return this->at(lor.index()); }

  /// Returns constant LOR geometry description for given LOR.
  const LORGeometry& operator[](const LOR& lor) const {
    return this->at(lor.index());
  }

  /// Append pixel info for given LOR.
  void push_back_pixel_info(const LOR& lor, const PixelInfo& pixel_info) {
    (*this)[lor].pixel_infos.push_back(pixel_info);
  }

  /// Append pixel weight for given LOR.
  void push_back_pixel(const LOR& lor, const Pixel& pixel, F weight) {
    auto& lor_geometry = (*this)[lor];
    auto center = grid.center_at(pixel.x, pixel.y);
    auto t = lor_geometry.segment.projection_scaled(center);
    PixelInfo pixel_info;
    pixel_info.pixel = pixel;
    pixel_info.t = t;
    pixel_info.weight = weight;
    push_back_pixel_info(lor, pixel_info);
  }

  void put_pixel(const LOR& lor, const Pixel& pixel, F weight) {
    auto& lor_geometry = (*this)[lor];
    lor_geometry.put_pixel(pixel, weight);
  }

  /// Sort all LOR geometry descriptions.
  ////
  /// This does not affect order of LOR gemetry descriptions here.
  void sort_all() {
    for (auto& lor_geometry : *this) {
      lor_geometry.sort();
    }
  }

  void sort_all_by_index() {
    for (auto& lor_geometry : *this) {
      lor_geometry.sort_by_index();
    }
  }

  /// Remove all pixel information from all LOR geometry descriptions.
  void erase_pixel_info() {
    for (auto& lor_geometry : *this) {
      lor_geometry.pixel_infos.resize(0);
    }
  }

  /// Return total number of pixel infos hold for all LORs.
  size_t n_pixel_infos() const {
    size_t total = 0;
    for (const auto& lor_geometry : *this) {
      total += lor_geometry.pixel_infos.size();
    }
    return total;
  }

  const S n_detectors;  ///< number of detectors
  const Grid grid;      ///< pixel grid description
};

}  // Barrel
}  // PET2D
