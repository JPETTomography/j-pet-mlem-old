#pragma once

#include <fstream>
#include <iomanip>

namespace util {

/// \c JSON file generator based on \c std::ofstream
///
/// Actually it does not currently nothing, but it is used to distinguish << >>
/// operators.
class json_ostream : public std::ofstream {
 public:
  /// Constructs new \a SVG file at given path with provided dimensions
  json_ostream(const std::string& fn  ///< Path to \a JSON file
               )
      : std::ofstream(fn) {}

  ~json_ostream() {}

  /// conditionally appends \c , delimiter or flips the argument flag
  json_ostream& delimiter(bool& next) {
    if (next) {
      *this << ", ";
    } else {
      next = true;
    }
    return *this;
  }
};
}  // util
