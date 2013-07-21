class png_writer {
 public:
  png_writer(std::string fn);

  template <typename byte_type = uint8_t>
  void write_header(unsigned int width, unsigned int height) {
    // only 8bpp & 16bpp are supported by PNG
    priv_write_header(width, height, sizeof(byte_type) > 1 ? 16 : 8);
  }

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
