#ifndef STRING_H
#define STRING_H

#include <string>
#include <list>

void split_string_on(std::string in,
                     std::string dels,
                     std::list<std::string>& parts);

#endif  // STRING_H
