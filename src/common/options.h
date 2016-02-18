#pragma once

#include "cmdline.h"

/// \cond PRIVATE
namespace Common {

void add_cuda_options(cmdline::parser& cl);
void add_openmp_options(cmdline::parser& cl);

template <typename T1, typename T2> struct Convert {
  static std::vector<T1> cast(const std::vector<T2>& vector) {
    std::vector<T1> ret;
    for (const auto& el : vector) {
      ret.push_back(el);
    }
    return ret;
  }
};

template <typename T> struct Convert<T, T> {
  static std::vector<T> cast(const std::vector<T>& vector) { return vector; }
};

}  // Common
/// \endcond
