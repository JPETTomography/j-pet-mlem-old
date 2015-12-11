#pragma once
#include <string>
#include <sstream>
#include <iomanip>

// \cond PRIVATE

// redefine help formatting for greater readibility
namespace cmdline {

class path : public std::string {
 public:
  path() : std::string() {}

  path(const char* s) : std::string(s) {}

  path(const std::string& str) : std::string(str) {}

  path(const std::string&& str) : std::string(str) {}

  path wo_ext() const { return substr(0, ext_pos(sep_pos())); }

  path wo_path() const {
    auto fn_sep = sep_pos();
    return substr(fn_sep != std::string::npos ? fn_sep + 1 : 0);
  }

  path ext() const {
    auto fn_ext = ext_pos(sep_pos());
    return substr(fn_ext != std::string::npos ? fn_ext : size(), size());
  }

  /// Add index suffix to file name with given zero padding
  path add_index(int index,             ///< index to append as suffix
                 int max_or_digits = 0  ///< maximum value
                                        ///  (or when negative number of digits)
                 ) const {
    // NOTE: This assumes index is <= max (not index < max).
    int n_digits = 1;
    if (max_or_digits > 0) {
      int max_for_digits = 10;
      while (max_for_digits <= max_or_digits) {
        ++n_digits;
        max_for_digits *= 10;
      }
    } else {
      n_digits = max_or_digits;
    }
    std::stringstream fn;
    fn << *this << "_" << std::setw(n_digits)  //
       << std::setfill('0')                    //
       << index << std::setw(0);               // 001
    return fn.str();
  }

  /// Return index suffix value and number of its digits (if present)
  int scan_index(int& n_digits) const {
    size_t start;
    for (start = length(); start > 0; --start) {
      if (!std::isdigit((*this)[start - 1]))
        break;
    }
    n_digits = length() - start;
    std::istringstream ss(substr(start));
    int index;
    ss >> index;
    return index;
  }

  /// Return index suffix value
  int scan_index() const {
    int n_digits;
    return scan_index(n_digits);
  }

 private:
  std::string::size_type sep_pos() const { return find_last_of("\\/"); }

  std::string::size_type ext_pos(const std::string::size_type fn_sep) const {
    auto fn_ext = find_last_of(".");
    if (fn_sep != std::string::npos && fn_sep > fn_ext)
      fn_ext = std::string::npos;
    return fn_ext;
  }
};

namespace detail {
template <> inline std::string readable_typename<int>() { return "int"; }
template <> inline std::string readable_typename<long>() { return "seed"; }
template <> inline std::string readable_typename<double>() { return "float"; }
template <> inline std::string readable_typename<path>() { return "file"; }

template <> inline std::string default_value<int>(int def) {
  if (def == 0)
    return "auto";
  return detail::lexical_cast<std::string>(def);
}
template <> inline std::string default_value<double>(double def) {
  if (def == 0)
    return "auto";
  return detail::lexical_cast<std::string>(def);
}
#if __GNUC__
template <> inline std::string default_value<ssize_t>(ssize_t def) {
  if (def < 0)
    return "all";
  return detail::lexical_cast<std::string>(def);
}
#endif

// teach cmdline handle k, m, g suffixes
namespace {
static inline long cast_to_long(const std::string& str) {
  const char* cstr = str.c_str();
  char* endstr;
  auto ret = std::strtol(cstr, &endstr, 10);
  if (endstr[0] && !endstr[1]) {
    switch (std::tolower(endstr[0])) {
      case 'k':
        ret *= 1000;
        break;
      case 'm':
        ret *= 1000000;
        break;
      case 'g':
        ret *= 1000000000;
        break;
      default:
        throw std::bad_cast();
    }
  } else if (endstr[0] && endstr[1])
    throw std::bad_cast();
  return ret;
}
}

template <> class lexical_cast_t<int, std::string, false> {
 public:
  static int cast(const std::string& arg) { return cast_to_long(arg); }
};

template <> class lexical_cast_t<long, std::string, false> {
 public:
  static int cast(const std::string& arg) { return cast_to_long(arg); }
};

}  // detail
}  // cmdline
// \endcond
