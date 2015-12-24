#pragma once

#include "util/bstream.h"
#include "util/read.h"

#include "2d/barrel/lor.h"
#include "2d/geometry/pixel_grid.h"
#include "2d/geometry/point.h"
#include "2d/geometry/line_segment.h"

namespace PET2D {
namespace Barrel {

/// Keeps geometry information about LOR
////
/// Keeps geometry information about LOR such as pixels belonging to the LOR
/// together with their description using PixelInfo.
template <typename FType, typename SType> struct LORGeometry {
  using F = FType;
  using S = SType;
  using Pixel = PET2D::Pixel<S>;
  using LOR = PET2D::Barrel::LOR<S>;
  using LineSegment = PET2D::LineSegment<F>;

  /// Information about given Pixel relative to the LOR
  struct PixelInfo {
    Pixel pixel;  ///< pixel itself
    F t;          ///< position of projected Pixel along the LOR
    F distance;   ///< distance of Pixel to LOR
    F fill;       ///< percentage of Pixel inside the LOR
    F weight;     ///< custom weight of the Pixel eg. probability
  };

  using PixelInfoList = std::vector<PixelInfo>;

  LOR lor;
  LineSegment segment;
  F width;
  PixelInfoList pixel_infos;

  LORGeometry() = default;

  LORGeometry(const LOR& lor, const LineSegment& segment, const F width)
      : lor(lor), segment(segment), width(width) {}

  /// Construct LOR info from stream.
  LORGeometry(std::istream& in)
      : lor(in), segment(in), width(util::read<F>(in)) {
    size_t n_pixels;
    in >> n_pixels;
    if (in && n_pixels) {
      pixel_infos.resize(n_pixels);
      in.read((char*)&pixel_infos[0], sizeof(PixelInfo) * n_pixels);
    }
  }

  /// Construct LOR info from binary stream.
  LORGeometry(util::ibstream& in) : lor(in), segment(in), width(in.read<F>()) {
    size_t n_pixels;
    in >> n_pixels;
    if (in && n_pixels) {
      pixel_infos.resize(n_pixels);
      in.read((char*)&pixel_infos[0], sizeof(PixelInfo) * n_pixels);
    }
  }

  /// Output LOR geometry to text stream.
  friend std::ostream& operator<<(std::ostream& out,
                                  const LORGeometry& lor_info) {
    out << lor_info.lor.first << ' ' << lor_info.lor.second << ' '
        << lor_info.segment << ' '  //
        << lor_info.width << ' '    //
        << 0;  // FIXME: unsupported serialization to text stream
    return out;
  }

  /// Output LOR geometry to binary stream.
  friend util::obstream& operator<<(util::obstream& out,
                                    const LORGeometry& lor_info) {
    out << lor_info.lor;
    out << lor_info.segment;
    out << lor_info.width;
    out << lor_info.pixel_infos.size();
    out << lor_info.pixel_infos;
    return out;
  }

  /// Sort pixel infos using position along LOR.
  void sort() {
    std::sort(pixel_infos.begin(),
              pixel_infos.end(),
              [](const PixelInfo& a, const PixelInfo& b) { return a.t < b.t; });
  }

  /// Sort pixel infos by pixel \c y then \c x.
  void sort_by_pixel() {
    std::sort(pixel_infos.begin(),
              pixel_infos.end(),
              [](const PixelInfo& a, const PixelInfo& b) {
                return a.pixel < b.pixel;
              });
  }

  /// Add weight for given pixel.
  ////
  /// Assumes that pixels in LOR are sorted by \c y then \c x.
  ///
  /// This is used to load external weights for geometry, eg. from existing
  /// Monte-Carlo matrix file.
  void add_weight_for_pixel(const Pixel& pixel, F weight) {
    PixelInfo reference_pixel_info;
    reference_pixel_info.pixel = pixel;
    auto pixel_info =
        std::lower_bound(pixel_infos.begin(),
                         pixel_infos.end(),
                         reference_pixel_info,
                         [](const PixelInfo& a, const PixelInfo& b) {
                           return a.pixel < b.pixel;
                         });

    // check if we found the right pixel
    if (pixel_info == pixel_infos.end() || pixel_info->pixel != pixel) {
#if THROW_ON_MISSING_PIXEL
      std::stringstream msg;
      msg << "pixel (" << pixel.x << "," << pixel.y << ") not found in LOR ";
      msg << lor.first << "-" << lor.second;
      for (const auto& pixel_info : pixel_infos) {
        msg << " (" << pixel_info.pixel.x << "," << pixel_info.pixel.y << ")";
      }
      throw msg.str();
#else
      return;
#endif
    }

#if NO_DUPLICATES
    if (pixel_info->weight > 0) {
      std::cerr << "duplicate entry for: lor " << lor.first << "-" << lor.second
                << " pixel " << pixel.x << "," << pixel.y <<  //
          " (" << pixel_info->pixel.x << "," << pixel_info->pixel.y << ")"
                << std::endl;
    }
#endif
    pixel_info->weight += weight;
  }

  /// Append given PixelInfo to the pixel info list.
  void push_back(const PixelInfo& pixel_info) {
    pixel_infos.push_back(pixel_info);
  }
};

}  // Barrel
}  // PET2D
