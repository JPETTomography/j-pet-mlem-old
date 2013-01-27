#pragma once

// @cond PRIVATE
// redefine help formatting for greater readibility
namespace cmdline {
class string : public std::string {
};
namespace detail {
template <> inline std::string readable_typename<ssize_t>() { return "index"; }
template <> inline std::string readable_typename<size_t>() { return "size"; }
template <> inline std::string readable_typename<double>() { return "float"; }
template <> inline std::string readable_typename<string>() { return "file"; }

template <> inline std::string default_value<double>(double def) {
  if (def == 0.)
    return "auto";
  return detail::lexical_cast<std::string>(def);
}
template <> inline std::string default_value<ssize_t>(ssize_t def) {
  if (def < 0)
    return "all";
  return detail::lexical_cast<std::string>(def);
}
}
}
// @end cond
