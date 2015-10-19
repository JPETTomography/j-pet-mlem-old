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

  path add_index(int index, int max) const {
    // NOTE: This assumes index is <= max (not index < max).
    int n_digits = 1;
    int max_for_digits = 10;
    while (max_for_digits <= max) {
      ++n_digits;
      max_for_digits *= 10;
    }
    std::stringstream fn;
    fn << *this << "_" << std::setw(n_digits)  //
       << std::setfill('0')                    //
       << index << std::setw(0);               // 001
    return fn.str();
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
template <> inline std::string readable_typename<int>() { return "size"; }
template <> inline std::string readable_typename<long>() { return "seed"; }
template <> inline std::string readable_typename<double>() { return "float"; }
template <> inline std::string readable_typename<path>() { return "file"; }

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
}

}  // cmdline
// \endcond
