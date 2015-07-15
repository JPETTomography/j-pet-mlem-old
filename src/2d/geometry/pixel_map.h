#pragma once

#include <cstdint>
#include <vector>
#include <ostream>

#include "util/png_writer.h"

namespace PET2D {

/// Rectangular map of pixels aka 2D image

/// Can be used to write pixels to PNG or generic stream
template <typename PixelType, typename ValueType>
class PixelMap : std::vector<ValueType> {
 public:
  using Pixel = PixelType;
  using S = typename Pixel::S;
  using Size = typename Pixel::Size;
  using Value = ValueType;
  using Base = std::vector<Value>;
  using iterator = typename Base::iterator;
  using const_iterator = typename Base::const_iterator;

  /// Creates new pixel map of given dimensions and default value
  PixelMap(S width, S height, ValueType value = 0)
      : Base(static_cast<Size>(width) * height, value),
        width(width),
        height(height) {}
  /// Copy constructor
  PixelMap(const PixelMap& other)
      : Base(other), width(other.width), height(other.height) {}
  /// Move constructor
  PixelMap(PixelMap&& other)
      : Base(other), width(other.width), height(other.height) {}

  PixelMap& operator=(const PixelMap& other) {
    if (other.width != width || other.height != height)
      throw("cannot assign pixel map of diffent size");
    Base::operator=(other);
    return *this;
  }
  PixelMap& operator=(PixelMap&& other) {
    if (other.width != width || other.height != height)
      throw("cannot assign pixel map of diffent size");
    Base::operator=(other);
    return *this;
  }

  Value& operator[](const Pixel& pixel) { return this->at(pixel.index(width)); }
  const Value& operator[](const Pixel& pixel) const {
    return this->at(pixel.index(width));
  }
  Value& operator[](Size index) { return this->at(index); }
  const Value& operator[](Size index) const { return this->at(index); }

  iterator begin() { return Base::begin(); }
  const_iterator begin() const { return Base::begin(); }
  iterator end() { return Base::end(); }
  const_iterator end() const { return Base::end(); }

  const S width;
  const S height;

  friend util::png_writer& operator<<(util::png_writer& png,
                                      const PixelMap& map) {
    png.write(map.width, map.height, map.data());
    return png;
  }

  friend std::ostream& operator<<(std::ostream& out, const PixelMap& map) {
    auto it = map.begin();
    auto end = map.end();
    auto width = map.width;
    for (int c = 1; it != end; ++it, ++c) {
      out << *it;
      if (c % width == 0) {
        out << "\n";
      } else {
        out << " ";
      }
    }
    return out;
  }
};

}  // PET2D
