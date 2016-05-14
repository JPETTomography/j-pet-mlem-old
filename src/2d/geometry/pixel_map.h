#pragma once

#if !__CUDACC__
#include <ostream>
#include "util/png_writer.h"
#include "util/nrrd_writer.h"
#include "util/bstream.h"
#endif

namespace PET2D {

/// Rectangular map of pixels aka 2D image
////
/// Can be used to write pixels to PNG or generic stream
template <typename PixelType, typename ValueType> class PixelMap {
 public:
  using Pixel = PixelType;
  using S = typename Pixel::S;
  using Size = typename Pixel::Size;
  using Value = ValueType;
  using iterator = Value*;
  using const_iterator = const Value*;

  /// Creates new pixel map of given dimensions initialized to zero
  PixelMap(S width, S height)
      : size(static_cast<Size>(width) * height),
        width(width),
        height(height),
        data(new Value[size]()),
        owning(true) {}

  /// Creates new pixel map of given dimensions and default value
  PixelMap(S width, S height, const Value& value)
      : size(static_cast<Size>(width) * height),
        width(width),
        height(height),
        data(new Value[size]),
        owning(true) {
    assign(value);
  }

  /// Creates new pixel map of given dimensions and existing data
  PixelMap(S width, S height, Value* data)
      : size(static_cast<Size>(width) * height),
        width(width),
        height(height),
        data(data),
        owning(false) {}

  /// Copy constructor
  PixelMap(const PixelMap& other)
      : size(other.size),
        width(other.width),
        height(other.height),
        data(other.owning ? new Value[size] : other.data),
        owning(other.owning) {
    if (owning) {
      memcpy(data, other.data, size * sizeof(Value));
    }
  }

  /// Move constructor
  PixelMap(PixelMap&& other)
      : size(other.size),
        width(other.width),
        height(other.height),
        data(other.data),
        owning(other.owning) {
    // make other object not to release its memory
    other.owning = false;
  }

  /// Destroys underlying data (only if it is managed)
  ~PixelMap() {
    if (owning && data) {
      delete[] data;
    }
  }

  Value& operator[](const Pixel& pixel) { return data[pixel.index(width)]; }
  const Value& operator[](const Pixel& pixel) const {
    return data[pixel.index(width)];
  }
  Value& operator[](Size index) { return data[index]; }
  const Value& operator[](Size index) const { return data[index]; }

  /// Increments atomically value
  void increment(const Pixel& pixel) {
#if _OPENMP && !_MSC_VER
    // FIXME: on MSVC it does not work
    __atomic_add_fetch(&data[pixel.index(width)], 1, __ATOMIC_SEQ_CST);
#else
    data[pixel.index(width)]++;
#endif
  }

  void assign(const Value& value) {
    for (auto& element : *this) {
      element = value;
    }
  }

  iterator begin() { return data; }
  const_iterator begin() const { return data; }
  iterator end() { return &data[size]; }
  const_iterator end() const { return &data[size]; }

  const Size size;
  const S width;
  const S height;
  Value* const data;

#if !__CUDACC__
  friend util::png_writer& operator<<(util::png_writer& png,
                                      const PixelMap& map) {
    png.write(map.width, map.height, map.data);
    return png;
  }

  friend util::nrrd_writer& operator<<(util::nrrd_writer& nrrd,
                                       const PixelMap& map) {
    size_t dimensions[2] = { (size_t)map.width, (size_t)map.height };
    nrrd.write_header<ValueType, 2>(dimensions);
    return nrrd;
  }

  friend std::istream& operator>>(std::istream& in, PixelMap& map) {
    for (auto voxel : map) {
      in >> voxel;
    }
    return in;
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

  friend util::ibstream& operator>>(util::ibstream& in, PixelMap& map) {
    in.read(map.data, map.size);
    return in;
  }

  friend util::obstream& operator<<(util::obstream& out, const PixelMap& map) {
    out.write(map.data, map.size);
    return out;
  }
#endif

 private:
  bool owning;  ///< are we owning the data pointer?
};

}  // PET2D
