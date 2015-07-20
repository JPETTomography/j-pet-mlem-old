#pragma once

#include <vector>
#include <ostream>

#include "util/png_writer.h"

namespace PET3D {

/// Cubical map of voxels aka 3D image

/// Can be used to write voxels to PNG or generic stream
template <typename VoxelType, typename ValueType>
class VoxelMap : std::vector<ValueType> {
 public:
  using Voxel = VoxelType;
  using S = typename Voxel::S;
  using Size = typename Voxel::Size;
  using Value = ValueType;
  using Base = std::vector<Value>;
  using iterator = typename Base::iterator;
  using const_iterator = typename Base::const_iterator;

  /// Creates new voxel map of given dimensions and default value
  VoxelMap(S width, S height, S depth, ValueType value = 0)
      : Base(static_cast<Size>(width) * height * depth, value),
        width(width),
        height(height) {}
  /// Copy constructor
  VoxelMap(const VoxelMap& other)
      : Base(other),
        width(other.width),
        height(other.height),
        depth(other.depth) {}
  /// Move constructor
  VoxelMap(VoxelMap&& other)
      : Base(other),
        width(other.width),
        height(other.height),
        depth(other.depth) {}

  VoxelMap& operator=(const VoxelMap& other) {
    if (other.width != width || other.height != height || other.depth != depth)
      throw("cannot assign voxel map of diffent size");
    Base::operator=(other);
    return *this;
  }
  VoxelMap& operator=(VoxelMap&& other) {
    if (other.width != width || other.height != height || other.depth != depth)
      throw("cannot assign voxel map of diffent size");
    Base::operator=(other);
    return *this;
  }

  Value& operator[](const Voxel& voxel) {
    return this->at(voxel.index(width, depth));
  }
  const Value& operator[](const Voxel& voxel) const {
    return this->at(voxel.index(width, depth));
  }
  Value& operator[](Size index) { return this->at(index); }
  const Value& operator[](Size index) const { return this->at(index); }

  void assign(const Value& v) { Base::assign(this->size(), v); }

  iterator begin() { return Base::begin(); }
  const_iterator begin() const { return Base::begin(); }
  iterator end() { return Base::end(); }
  const_iterator end() const { return Base::end(); }

  const S width;
  const S height;
  const S depth;

  friend util::png_writer& operator<<(util::png_writer& png,
                                      const VoxelMap& map) {
    png.write(map.width, map.height * map.depth, map.data());
    return png;
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
};

}  // PET2D
