#include <string>
#include <ostream>

#ifdef HAVE_LIBPNG
#include <png.h>
struct png_writer_private {
  FILE *fp;
  png_structp png_ptr;
  png_infop info_ptr;
};
#endif

#include "png_writer.h"

png_writer::png_writer(std::string fn) {
#ifdef HAVE_LIBPNG
  priv = new png_writer_private;
  priv->fp       = nullptr;
  priv->png_ptr  = nullptr;
  priv->info_ptr = nullptr;

  if (!( priv->png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL) )) {
    throw(std::string("cannot create png version"));
  }

  if (!( priv->info_ptr = png_create_info_struct(priv->png_ptr) )) {
    throw(std::string("cannot create png info"));
  }

  if (setjmp(png_jmpbuf(priv->png_ptr))) {
    throw(std::string("cannot hook png exception"));
  }

  if (!( priv->fp = fopen(fn.c_str(), "wb") )) {
    throw(std::string("cannot create output file"));
  }

  png_init_io(priv->png_ptr, priv->fp);
#endif
}

void png_writer::priv_write_header(size_t width, size_t height, size_t bpp) {
#ifdef HAVE_LIBPNG
  png_set_IHDR(priv->png_ptr, priv->info_ptr, width, height, bpp,
        PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE,
        PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
  png_write_info(priv->png_ptr, priv->info_ptr);
#endif
}

void png_writer::priv_write_row(unsigned char *row) {
#ifdef HAVE_LIBPNG
  png_write_row(priv->png_ptr, row);
#endif
}

png_writer::~png_writer() {
#ifdef HAVE_LIBPNG
  if (priv->fp) {
    png_write_end(priv->png_ptr, NULL);
    fclose(priv->fp);
  }
  if (priv->info_ptr) png_free_data(priv->png_ptr, priv->info_ptr, PNG_FREE_ALL, -1);
  if (priv->png_ptr) png_destroy_write_struct(&priv->png_ptr, (png_infopp)NULL);
  if (priv) delete priv;
#endif
}
