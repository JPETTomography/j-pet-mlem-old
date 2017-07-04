#include <string>
#include <iostream>

#include "util/string.h"

void split_string_on(std::string in,
                     std::string dels,
                     std::list<std::string>& parts) {
  std::string::size_type start = 0;
  while (true) {
    start = in.find_first_not_of(dels, start);
    if (start == std::string::npos)
      return;

    auto end = in.find_first_of(dels, start);
    if (end == std::string::npos) {
      parts.push_back(in.substr(start));
      return;
    } else {
      parts.push_back(in.substr(start, end - start));
      start = end;
    }
  }
}
