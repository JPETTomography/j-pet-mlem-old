#pragma once

#if !__CUDACC__
#include <ostream>
#include "util/png_writer.h"
#include "util/nrrd_writer.h"
#include "util/bstream.h"
#endif

namespace PET3D {

/// Cubical map of voxels aka 3D image
////
/// Can be used to write voxels to PNG or generic stream
template <typename VoxelType, typename ValueType> class VoxelMap {
 public:
  using Voxel = VoxelType;
  using S = typename Voxel::S;
  using Size = typename Voxel::Size;
  using Value = ValueType;
  using iterator = Value*;
  using const_iterator = const Value*;

  /// Creates new voxel map of given dimensions initialized to zero
  VoxelMap(S width, S height, S depth)
      : size(static_cast<Size>(width) * height * depth),
        width(width),
        height(height),
        depth(depth),
        data(new Value[size]()),
        owning(true) {}

  /// Creates new voxel map of given dimensions and default value
  VoxelMap(S width, S height, S depth, const Value& value)
      : size(static_cast<Size>(width) * height * depth),
        width(width),
        height(height),
        depth(depth),
        data(new Value[size]),
        owning(true) {
    assign(value);
  }

  /// Creates new voxel map of given dimensions and existing data
  VoxelMap(S width, S height, S depth, Value* data)
      : size(static_cast<Size>(width) * height * depth),
        width(width),
        height(height),
        depth(depth),
        data(data),
        owning(false) {}

  /// Copy constructor
  VoxelMap(const VoxelMap& other)
      : size(other.size),
        width(other.width),
        height(other.height),
        depth(other.depth),
        data(other.owning ? new Value[size] : other.data),
        owning(other.owning) {
    if (owning) {
      memcpy(data, other.data, size * sizeof(Value));
    }
  }

  /// Move constructor
  VoxelMap(VoxelMap&& other)
      : size(other.size),
        width(other.width),
        height(other.height),
        depth(other.depth),
        data(other.data),
        owning(other.owning) {
    // make other object not to release its memory
    other.owning = false;
  }

  /// Destroys underlying data (only if it is managed)
  ~VoxelMap() {
    if (owning && data) {
      delete[] data;
    }
  }

  Value& operator[](const Voxel& voxel) {
    return data[voxel.index(width, height)];
  }
  const Value& operator[](const Voxel& voxel) const {
    return data[voxel.index(width, height)];
  }
  Value& operator[](Size index) { return data[index]; }
  const Value& operator[](Size index) const { return data[index]; }

  void assign(const Value& value) {
    for (auto& element : *this) {
      element = value;
    }
  }

  /// Extract from other voxel map at source origin relative to source center
  void copy(const VoxelMap& source, const Voxel origin) {
    Voxel translation(source.width, source.height, source.depth);
    translation -= Voxel(width, height, depth);
    translation.x /= 2, translation.y /= 2, translation.z /= 2;
    translation += origin;
    for (S z = 0; z < size; ++z) {
      for (S y = 0; y < size; ++y) {
        for (S x = 0; x < size; ++x) {
          Voxel voxel(x, y, z);
          (*this)[voxel] = source[voxel + translation];
        }
      }
    }
  }

  iterator begin() { return data; }
  const_iterator begin() const { return data; }
  iterator end() { return &data[size]; }
  const_iterator end() const { return &data[size]; }

  const Size size;
  const S width;
  const S height;
  const S depth;
  Value* const data;

#if !__CUDACC__
  friend util::png_writer& operator<<(util::png_writer& png,
                                      const VoxelMap& map) {
    png.write(map.width, map.height * map.depth, map.data);
    return png;
  }

  friend util::nrrd_writer& operator<<(util::nrrd_writer& nrrd,
                                       const VoxelMap& map) {
    size_t dimensions[3] = {
      (size_t)map.width, (size_t)map.height, (size_t)map.depth
    };
    nrrd.write_header<ValueType, 3>(dimensions);
    return nrrd;
  }

  friend std::istream& operator>>(std::istream& in, VoxelMap& map) {
    for (auto voxel : map) {
      in >> voxel;
    }
    return in;
  }

  friend std::ostream& operator<<(std::ostream& out, const VoxelMap& map) {
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

  friend util::ibstream& operator>>(util::ibstream& in, VoxelMap& map) {
    in.read(map.data, map.size);
    return in;
  }

  friend util::obstream& operator<<(util::obstream& out, const VoxelMap& map) {
    out.write(map.data, map.size);
    return out;
  }
#endif

 private:
  bool owning;  ///< are we owning the data pointer?
};

}  // PET2D
