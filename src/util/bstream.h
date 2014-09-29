#pragma once

#include <fstream>
#include <vector>

class ibstream : public std::ifstream {
 public:
  ibstream(std::string fn, std::ios_base::openmode mode = std::ios_base::in)
      : std::ifstream(fn, mode | std::ios_base::in) {}

  template <typename T> ibstream& operator>>(T& v) {
    read(reinterpret_cast<char*>(&v), sizeof(v));
    return *this;
  }
};

class obstream : public std::ofstream {
 public:
  obstream(std::string fn, std::ios_base::openmode mode = std::ios_base::out)
      : std::ofstream(fn, mode | std::ios_base::out) {}

  template <typename T> obstream& operator<<(const T v) {
    write(reinterpret_cast<const char*>(&v), sizeof(v));
    return *this;
  }

  template <typename T> obstream& operator<<(const std::vector<T> vector) {
    for (auto&& v : vector) {
      *this << v;
    }
    return *this;
  }
};
