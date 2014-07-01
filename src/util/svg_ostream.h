#pragma once

#include <fstream>

template <typename FType = double> class svg_ostream : public std::ofstream {
 public:
  typedef FType F;

  svg_ostream(const std::string& fn,
              F x_max,
              F y_max,
              F image_width,
              F image_height)
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
    *this << "    polygon, .detector {"
          << " fill: #f99;"
          << " stroke: red;"
          << " }" << std::endl;
    *this << "    circle  { fill: white;"
          << " stroke: green;"
          << " }" << std::endl;
    *this << "    .circle_detector { fill: #ddd;"
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

  svg_ostream& link_image(std::string fn, F x, F y, F width, F height) {
    *this << "<image"
          << " xlink:href=\"" << fn << "\""
          << " x=\"" << x << "\""
          << " y=\"" << y << "\""
          << " height=\"" << width << "\""
          << " width=\"" << height << "\"/>" << std::endl;
    return *this;
  }
};
