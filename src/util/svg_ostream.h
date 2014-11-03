#pragma once

#include <fstream>

namespace util {

/// \c SVG file generator based on \c std::ofstream
template <typename FType = double> class svg_ostream : public std::ofstream {
 public:
  using F = FType;

  /// Constructs new \a SVG file at given path with provided dimensions
  svg_ostream(const std::string& fn,  ///< Path to \a SVG file
              F x_max,                ///< Maximum \c x value
              F y_max,                ///< Maximum \c y value
              F image_width,          ///< Image (canvas) width
              F image_height          ///< Image (canvas) height
              )
      : std::ofstream(fn) {
    auto x_translate = x_max;
    auto y_translate = y_max;
    auto scale = std::min(static_cast<F>(image_width - 4) / x_max / 2,
                          static_cast<F>(image_height - 4) / y_max / 2);
    auto stroke = 1. / scale;
    *this << "<?xml version=\"1.0\" standalone=\"no\"?>" << std::endl;
    *this << "<svg"
          << " width=\"" << image_width << "\""
          << " height=\"" << image_height << "\""
          << " version=\"1.1\""
          << " xmlns=\"http://www.w3.org/2000/svg\""
          << " xmlns:xlink=\"http://www.w3.org/1999/xlink\">" << std::endl;
    *this << "<defs>" << std::endl;
    *this << "  <style type=\"text/css\">"
          << "<![CDATA[" << std::endl;
#if BLACK_BACKGROUND
    *this << "    svg     { background: black; }" << std::endl;
#endif
    *this << "    polygon, circle {"
          << " stroke-width: 0.5;"
          << " vector-effect: non-scaling-stroke;"
          << " }" << std::endl;
    *this << "    @media print { polygon, circle {"
          << " stroke-width: " << stroke << ";"
          << " vector-effect: none;"
          << " } }" << std::endl;
    *this << "    polygon, #scintillators circle {"
          << " fill: #f99;"
          << " stroke: red;"
          << " }" << std::endl;
    *this << "    circle  { fill: white;"
          << " stroke: green;"
          << " }" << std::endl;
    *this << "    #photomultipiers circle { fill: #ddd;"
          << " stroke: #999;"
          << " }" << std::endl;
    *this << "  ]]>"
          << "</style>" << std::endl;
    *this << "</defs>" << std::endl;
    *this << "<g transform=\"translate(2, 2)\">" << std::endl;
    *this << "<g transform=\"scale(" << scale << ',' << scale << ")\">"
          << std::endl;
    *this << "<g transform=\"translate(" << x_translate << "," << y_translate
          << ")\">" << std::endl;
  }

  ~svg_ostream() {
    *this << "</g>" << std::endl;
    *this << "</g>" << std::endl;
    *this << "</g>" << std::endl;
    *this << "</svg>" << std::endl;
  }

  /// Embeds image at given path at provided position and dimensions
  svg_ostream& link_image(std::string fn,  ///< Path to embedded image
                          F x,             ///< Image \c x position
                          F y,             ///< Image \c y position
                          F width,         ///< Image width
                          F height         ///< Image height
                          ) {
    *this << "<image"
          << " xlink:href=\"" << fn << "\""
          << " x=\"" << x << "\""
          << " y=\"" << y << "\""
          << " height=\"" << width << "\""
          << " width=\"" << height << "\"/>" << std::endl;
    return *this;
  }
};
}  // util
