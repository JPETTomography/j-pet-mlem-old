#pragma once

#include <fstream>
#include <sstream>

namespace util {

/// Writes image header to \a NRRD file
////
/// \see http://teem.sourceforge.net/nrrd/
class nrrd_writer : std::ofstream {
 public:
  /// Constructs new \a NRRD writer writing to given path
  nrrd_writer(std::string fn, std::string raw_fn, bool textual = false)
      : std::ofstream(fn), raw_fn(raw_fn), textual(textual) {}

  /// Writes \a NRRD file header defining image dimensions and type
  template <typename value_type = uint8_t, size_t n_dimensions>
  void write_header(const size_t dimensions[n_dimensions]) {
    *this << "NRRD0001" << std::endl
          << "type: " << out_type<value_type>() << std::endl
          << "endian: little" << std::endl
          << "dimension: " << n_dimensions << std::endl
          << "sizes: " << out_dimensions(dimensions, n_dimensions) << std::endl
          << "data file: ./" << raw_fn << std::endl
          << "encoding: " << (textual ? "ascii" : "raw") << std::endl;
  }

 private:
  std::string out_dimensions(const size_t* dimensions, size_t n_dimensions) {
    std::stringstream ss;
    for (size_t i = 0; i < n_dimensions; ++i) {
      if (i > 0) {
        ss << " ";
      }
      ss << dimensions[i];
    }
    return ss.str();
  }

  template <typename value_type> const char* out_type() {
    throw("this type is not supported by NRRD");
    return nullptr;
  }

  const std::string raw_fn;
  bool textual;
};

template <> inline const char* nrrd_writer::out_type<float>() {
  return "float";
}
template <> inline const char* nrrd_writer::out_type<double>() {
  return "double";
}

}  // util
