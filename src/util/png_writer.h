#pragma once

#include <cstdint>

namespace util {

/// Writes custom image data to \a PNG file
class png_writer {
 public:
  /// Constructs new \a PNG writer writing to given path
  png_writer(std::string fn);

  /// Writes \a PNG file header defining image dimensions
  template <typename byte_type = uint8_t>
  void write_header(unsigned int width, unsigned int height) {
    // only 8bpp & 16bpp are supported by PNG
    priv_write_header(width, height, sizeof(byte_type) > 1 ? 16 : 8);
  }

  /// Writes image row

  /// Provided array should match dimensions provided by \ref write_header
  template <typename byte_type = uint8_t> void write_row(byte_type* row) {
    priv_write_row(reinterpret_cast<unsigned char*>(row));
  }

  ~png_writer();

 private:
  void priv_write_header(unsigned int width,
                         unsigned int height,
                         unsigned int bpp);
  void priv_write_row(unsigned char* row);
  struct png_writer_private* priv;
};
}  // util
