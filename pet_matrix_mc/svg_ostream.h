#pragma once

#include <fstream>

template <typename F = double>
class svg_ostream : public std::ofstream {
public:
  svg_ostream(const std::string fn, F x_max, F y_max, F image_width, F image_height)
  : std::ofstream(fn)
  {
    *this << "<?xml version=\"1.0\" standalone=\"no\"?>" << std::endl;
    *this << "<svg width=\"" << image_width
           << "\" height=\"" << image_height
           << "\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">"
          << std::endl;
    auto x_translate = x_max;
    auto y_translate = y_max;
    auto scale = std::min((image_width-4.) / x_max / 2., (image_height-4.) / y_max / 2.);
    *this << "<g transform=\"translate(2, 2)\">" << std::endl;
    *this << "<g transform=\"scale(" << scale << ")\">" << std::endl;
    *this << "<g transform=\"translate(" << x_translate
                                  << "," << y_translate
                                  << ")\">" << std::endl;
  }
  ~svg_ostream() {
    *this << "</g>" << std::endl;
    *this << "</g>" << std::endl;
    *this << "</g>" << std::endl;
    *this << "</svg>" << std::endl;
  }
};
