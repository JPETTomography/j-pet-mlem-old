#include "test.h"

#include "util/string.h"

TEST("split_string") {
  using string_list = std::list<std::string>;

  SECTION("one part") {
    std::string input = "identifier_1";
    string_list parts;
    split_string_on(input, " \t\n", parts);

    CHECK(parts.size() == 1);
    auto it = parts.begin();
    CHECK("identifier_1" == *it);
  }

  SECTION("one part leading del") {
    std::string input = " identifier_1";
    string_list parts;
    split_string_on(input, " \t\n", parts);

    CHECK(parts.size() == 1);
    auto it = parts.begin();
    CHECK("identifier_1" == *it);
  }

  SECTION("empty") {
    std::string input = " \t ";
    string_list parts;
    split_string_on(input, " \t\n", parts);

    CHECK(parts.size() == 0);
  }

  SECTION("one part leading del trailing del") {
    std::string input = " identifier_1\t\n";
    string_list parts;
    split_string_on(input, " \t\n", parts);

    CHECK(parts.size() == 1);
    auto it = parts.begin();
    CHECK("identifier_1" == *it);
  }

  SECTION("two parts leading del trailing del") {
    std::string input = " identifier_1 identifier_2\t\n";
    string_list parts;
    split_string_on(input, " \t\n", parts);

    CHECK(parts.size() == 2);
    auto it = parts.begin();
    CHECK("identifier_1" == *it);
    it++;
    CHECK("identifier_2" == *it);
  }

  SECTION("two parts") {
    std::string input = "identifier_1 identifier_2";
    string_list parts;
    split_string_on(input, " \t\n", parts);

    CHECK(parts.size() == 2);
    auto it = parts.begin();
    CHECK("identifier_1" == *it);
    it++;
    CHECK("identifier_2" == *it);
  }
}
