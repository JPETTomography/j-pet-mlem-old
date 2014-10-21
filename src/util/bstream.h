#pragma once

#include <fstream>
#include <vector>

namespace util {

/// Binary input stream based on \c std::ifstream

/// Reads any type using its native binary representation.
/// \note
/// This class should not be used to read pointers such as \c char*.
class ibstream : public std::ifstream {
 public:
  ibstream(std::string fn, std::ios_base::openmode mode = std::ios_base::in)
      : std::ifstream(fn, mode | std::ios_base::in) {}

  template <typename T> ibstream& operator>>(T& v) {
    std::ifstream::read(reinterpret_cast<char*>(&v), sizeof(v));
    return *this;
  }

  template <typename T> ibstream& read(T* ptr, size_t size) {
    std::ifstream::read(reinterpret_cast<char*>(ptr), sizeof(*ptr) * size);
    return *this;
  }
};

/// Binary output stream based on \c std::ifstream

/// Writes any type using its native binary representation.
/// \note
/// This class should not be used to read pointers such as \c char*.
class obstream : public std::ofstream {
 public:
  obstream(std::string fn, std::ios_base::openmode mode = std::ios_base::out)
      : std::ofstream(fn, mode | std::ios_base::out) {}

  template <typename T> obstream& operator<<(const T v) {
    std::ofstream::write(reinterpret_cast<const char*>(&v), sizeof(v));
    return *this;
  }

  template <typename T> obstream& operator<<(const std::vector<T>& vector) {
    for (auto&& v : vector) {
      *this << v;
    }
    return *this;
  }

  template <typename T> obstream& write(const T* ptr, size_t size) {
    std::ofstream::write(reinterpret_cast<const char*>(ptr),
                         sizeof(*ptr) * size);
    return *this;
  }
};
}  // util
