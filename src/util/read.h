#pragma once

#if !__CUDACC__
#include <istream>
namespace util {
template <typename Type> Type read(std::istream& in) {
  Type value;
  in >> value;
  return value;
}
}  // util
#endif
