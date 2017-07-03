#ifndef GATE_PARSER_H
#define GATE_PARSER_H

#include <stdio.h>

#include <string>
#include <list>

namespace Gate {

class Parser {
  class CommandChain {
   public:
    using string_list = std::list<std::string>;
    string_list::const_iterator commands() { return commands_.begin(); }

    friend class Parser;

   private:
    std::list<std::string> commands_;
    std::list<std::string> arguments_;
  };

 public:
  using CommandChain = Parser::CommandChain;

  CommandChain parse(std::string line) {
    CommandChain chain;
    auto start = line.find_first_of("/", 0);
    if (start == std::string::npos) {
      fprintf(stderr, "`%s' is not a gate command ", line.c_str());
    } else {
      start++;
      bool another_command = true;
      while (another_command) {
        auto next_backslash = line.find_first_of("/ \t\n", start);
        if (next_backslash == std::string::npos) {
          fprintf(stderr, "`%s' is not a gate command ", line.c_str());
          another_command = false;
        } else {
          auto command = line.substr(start, next_backslash - start);
          chain.commands_.push_back(command);
          if (line[next_backslash] == '/') {
            start = next_backslash + 1;
          } else {
            another_command = false;
          }
        }
      }
    }
    return chain;
  }
};
}

#endif  // GATE_PARSER_H
