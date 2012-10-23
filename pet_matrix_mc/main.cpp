// PET System Matrix Calculator
// Authors:
//   Adam Strzelecki <adam.strzlecki@uj.edu.pl>
//   Piotr Bialas <piotr.bialas@uj.edu.pl>
//
// Using Monte Carlo method and square detectorilators

#include <random>
#include <iostream>

#include <cmdline.h>
#include <png.h>

#include "detector_ring.h"

template <typename F> F deg(F rad) { return rad * 180/M_PI; }
template <typename F> F rad(F deg) { return deg * M_PI/180; }

namespace cmdline {
  namespace detail {
    template <> inline std::string readable_typename<ssize_t>() { return "index"; }
    template <> inline std::string readable_typename<size_t>()  { return "size"; }
    template <> inline std::string readable_typename<double>()  { return "float"; }

    template <>
    inline std::string default_value<double>(double def) {
      if (def == 0.) return "auto";
      return detail::lexical_cast<std::string>(def);
    }
    template <>
    inline std::string default_value<ssize_t>(ssize_t def) {
      if (def < 0) return "all";
      return detail::lexical_cast<std::string>(def);
    }
  }
}

template <typename F = double>
class always_accept {
public:
  always_accept() {}
  bool operator () (F) { return true; }
};

int main(int argc, char *argv[]) {

  cmdline::parser cl;

  cl.add<size_t>     ("n-pixels",    'n', "number of pixels in one dimension", false, 256);
  cl.add<size_t>     ("n-detectors", 'd', "number of ring detectors",          false, 64);
  cl.add<size_t>     ("n-emissions", 'e', "emissions per pixel",               false, 1);
  cl.add<double>     ("radious",     'r', "inner detector ring radious",       false);
  cl.add<double>     ("s-pixel",     'p', "pixel size",                        false);
  cl.add<double>     ("w-detector",  'w', "detector width",                    false);
  cl.add<double>     ("h-detector",  'h', "detector height",                   false);
  cl.add             ("stats",       's', "show stats");
  cl.add             ("wait",       '\0', "wait before exit");
  cl.add<ssize_t>    ("lor",         'l', "select lor to output to a file",    false, -1);
  cl.add<std::string>("output",      'o', "output a file",                     false);

  cl.parse_check(argc, argv);

  auto n_pixels    = cl.get<size_t>("n-pixels");
  auto n_detectors = cl.get<size_t>("n-detectors");
  auto n_emissions = cl.get<size_t>("n-emissions");
  auto radious     = cl.get<double>("radious");
  auto s_pixel     = cl.get<double>("s-pixel");
  auto w_detector  = cl.get<double>("w-detector");
  auto h_detector  = cl.get<double>("h-detector");

  // automatic pixel size
  if (radious == 0.) {
    if (cl.get<double>("s-pixel") == 0.) {
      radious = sqrt(2.);
    } else {
      radious = sqrt(s_pixel * n_pixels);
    }
    std::cerr << "--radious=" << radious << std::endl;
  }

  // automatic radious
  if (s_pixel == 0.) {
    if (cl.get<double>("radious") == 0.) {
      s_pixel = 2./n_pixels;
    } else {
      s_pixel = radious*radious / n_pixels;
    }
    std::cerr << "--s-pixel=" << s_pixel << std::endl;
  }

  // automatic detector size
  if (w_detector == 0.) {
    w_detector = 2 * M_PI * .9 * radious / n_detectors;
    std::cerr << "--w-detector=" << w_detector << std::endl;
  }
  if (h_detector == 0.) {
    h_detector = w_detector;
    std::cerr << "--h-detector=" << h_detector << std::endl;
  }

  std::random_device rd;
  std::mt19937 gen(rd());
  detector_ring<> dr(n_detectors, n_pixels, s_pixel, radious, w_detector, h_detector);
  dr.matrix_mc(gen, always_accept<>(), n_emissions, true, true);
  
  auto pixel_max = 0;
  auto pixel_min = std::numeric_limits<decltype(pixel_max)>::max();
  auto lor       = cl.get<ssize_t>("lor");

  if (cl.exist("stats") || cl.exist("output"))
    for (auto y = 0; y < n_pixels; ++y) {
      for (auto x = 0; x < n_pixels; ++x) {
        auto hits = lor > 0 ? dr.matrix(lor, x, y) : dr.hits(x, y);
        pixel_max = std::max(pixel_max, hits);
        pixel_min = std::min(pixel_min, hits);
      }
    }

  if (cl.exist("stats")) {
    std::cerr
      << "Non zero LORs: "
      << dr.non_zero_lors()
      << '/'
      << dr.lors()
      << std::endl;
    std::cerr
      << "Min hits: "
      << pixel_min
      << std::endl;
    std::cerr
      << "Max hits: "
      << pixel_max
      << std::endl;
  }

#ifdef HAVE_LIBPNG
  auto output = cl.get<std::string>("output");
  if (output.size()) {
    FILE *fp = nullptr;
    png_structp png_ptr = nullptr;
    png_infop info_ptr = nullptr;

    if (!( png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL) )) {
      throw(std::string("cannot create png version"));
      goto cleanup;
    }

    if (!( info_ptr = png_create_info_struct(png_ptr) )) {
      throw(std::string("cannot create png info"));
      goto cleanup;
    }

    if (setjmp(png_jmpbuf(png_ptr))) {
      throw(std::string("cannot hook png exception"));
      goto cleanup;
    }

    if (!( fp = fopen(output.c_str(), "wb") )) {
      throw(std::string("cannot create output file"));
      goto cleanup;
    }

    png_init_io(png_ptr, fp);

    // 16 bit gray
    png_set_IHDR(png_ptr, info_ptr, n_pixels, n_pixels,
          16, PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE,
          PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    png_write_info(png_ptr, info_ptr);

    for (auto y = 0; y < n_pixels; ++y) {
      uint16_t row[n_pixels];
      for (auto x = 0; x < n_pixels; ++x) {
        row[x] = static_cast<double>(lor > 0 ? dr.matrix(lor, x, y) : dr.hits(x, y))
          * std::numeric_limits<uint16_t>::max() / pixel_max;
      }
      png_write_row(png_ptr, reinterpret_cast<png_bytep>(row));
    }
    png_write_end(png_ptr, NULL);

  cleanup:
    if (fp) fclose(fp);
    if (info_ptr) png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
    if (png_ptr) png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
  }
#endif

  if (cl.exist("wait")) {
    std::cerr << "Press Enter." << std::endl;
    while(getc(stdin) != '\n');
  }

  return 0;
}
