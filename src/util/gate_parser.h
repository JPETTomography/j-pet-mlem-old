#ifndef GATE_PARSER_H
#define GATE_PARSER_H

#include <stdio.h>

#include <string>
#include <list>

namespace Gate {

class Parser {
  class CommandChain {
    CommandChain() : is_valid_(true) {}

   public:
    using string_list = std::list<std::string>;
    string_list::const_iterator commands() const { return commands_.begin(); }
    string_list::const_iterator arguments() const { return arguments_.begin(); }
    bool is_valid() const { return is_valid_; }

    friend class Parser;

   private:
    bool is_valid_;
    std::list<std::string> commands_;
    std::list<std::string> arguments_;
  };

  void parse_command(std::string input, CommandChain& chain) {
    auto start = input.find_first_of("/", 0);

    if (start == std::string::npos) {
      fprintf(stderr, "`%s' is not a gate command ", input.c_str());
      chain.is_valid_ = false;
      return;
    }

    start++;
    bool another_command = true;
    while (1) {
      auto next_backslash = input.find_first_of("/ \t\n", start);
      if (next_backslash == std::string::npos) {
        fprintf(stderr, "`%s' is not a gate command ", input.c_str());
        chain.is_valid_ = false;
        another_command = false;
        return;
      }
      auto command = input.substr(start, next_backslash - start);
      chain.commands_.push_back(command);
      if (input[next_backslash] == '/')
        start = next_backslash + 1;
      else
        break;
    }
  }

 public:
  using CommandChain = Parser::CommandChain;

  CommandChain parse(std::string line) {
    CommandChain chain;

    parse_command(line, chain);

    return chain;
  }
};
}

#endif  // GATE_PARSER_H
