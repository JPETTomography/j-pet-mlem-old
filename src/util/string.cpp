#include <string>
#include <iostream>

#include "util/string.h"

void split_string_on(std::string in,
                     std::string dels,
                     std::string::size_type pos,
                     std::list<std::string>& parts) {
  std::string::size_type start = pos;
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

void split_string_on(std::string in,
                     std::string dels,
                     std::list<std::string>& parts) {
  split_string_on(in, dels, 0, parts);
}
