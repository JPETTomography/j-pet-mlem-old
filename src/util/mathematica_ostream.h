#pragma once

#include <fstream>

namespace util {

/// Mathematica \c .m file generator based on \c std::ofstream

/// Actually it does not currently nothing, but it is used to distinguish << >>
/// operators.
class mathematica_ostream : public std::ofstream {
 public:
  /// Constructs new \a SVG file at given path with provided dimensions
  mathematica_ostream(const std::string& fn  ///< Path to \a .m file
                      )
      : std::ofstream(fn) {}

  ~mathematica_ostream() {}

  /// conditionally appends \c , delimiter or flips the argument flag
  mathematica_ostream& delimiter(bool& next) {
    if (next) {
      *this << ", ";
    } else {
      next = true;
    }
    return *this;
  }
};
}  // util
