#pragma once

#include <fcntl.h>
#include <unistd.h>
#include <cstring>

#include "delegate.h"

namespace util {

/// Fast text file parser, about 5x faster than reading with STL
////
/// Inspired and rewritten from:
/// \see https://github.com/alexflint/fscanf-speed
struct text_parser {
  /// Create parser for given text
  text_parser(const char* str) : str(str) {}

  /// Read given value from text file
  template <typename Type> text_parser& operator>>(Type&);

  /// Read quickly lines from file
  ////
  /// It reads 20 million lines in about 0.5 sec.
  static void read_lines(char const* fn,
                         delegate<void(const char* line)> got_line) {
    static const auto BUFFER_SIZE = 16 * 1024;
    int fd = ::open(fn, O_RDONLY);
    if (fd == -1)
      throw(std::string("cannot open input file: ") + fn);

    char buf[BUFFER_SIZE];
    size_t prev_size = 0;
    while (size_t bytes_read =
               ::read(fd, buf + prev_size, sizeof(buf) - prev_size)) {
      if (bytes_read == (size_t)-1)
        throw("read failed");
      if (!bytes_read)
        break;
      const auto buf_end = buf + prev_size + bytes_read;
      char *line, *line_end;
      for (line = buf; (line_end = (char*)memchr(line, '\n', buf_end - line));
           line = line_end + 1) {
        *line_end = 0;
        got_line(line);
      }
      prev_size = buf_end - line;
      if (prev_size > 0) {
        std::memmove(buf, line, prev_size);
      }
    }
  }

 private:
  const char* str;

  void skip_space() {
    while (*str == ' ')
      ++str;
  }

  static constexpr int kMaxExponent = 300;

  template <typename FType> static FType pow10_positive(int e) {
    static FType cache[kMaxExponent] = { (FType)1 };
    auto res = cache[e];
    if (!res) {
      res = (FType)10 * pow10_negative<FType>(e - 1);
      cache[e] = res;
    }
    return res;
  }

  template <typename FType> static FType pow10_negative(int e) {
    static FType cache[kMaxExponent] = { (FType)1 };
    auto res = cache[e];
    if (!res) {
      res = (FType)0.1 * pow10_negative<FType>(e - 1);
      cache[e] = res;
    }
    return res;
  }

  static int parse_integer(const char*& p, int& n_digits) {
    n_digits = 0;
    int sum = 0;
    while (*p >= '0' && *p <= '9') {
      sum = sum * 10 + (*p - '0');
      ++p;
      ++n_digits;
    }
    return sum;
  }

  static int parse_integer(const char*& p) {
    int n_digits;
    return parse_integer(p, n_digits);
  }

  static int parse_sign(const char*& p) {
    if (*p == '-') {
      ++p;
      return -1;
    } else if (*p == '+') {
      ++p;
      return 1;
    }
    return 1;
  }

  static int parse_signed_integer(const char*& p, int& n_digits) {
    auto sign = parse_sign(p);
    return sign * parse_integer(p, n_digits);
  }

  static int parse_signed_integer(const char*& p) {
    int n_digits;
    return parse_signed_integer(p, n_digits);
  }

  template <typename FType> static FType parse_fp(const char*& p) {
    // Consume whole part
    int n_digits_whole;
    int sign = parse_sign(p);
    int whole_part = parse_integer(p, n_digits_whole);
    // Consume fractional part
    FType val = whole_part;
    if (*p == '.') {
      ++p;
      int ndigits_fractional = 0;
      int fractional_part = parse_integer(p, ndigits_fractional);
      if (n_digits_whole == 0 && ndigits_fractional == 0) {
        throw("period with no digits either before or after");
      } else if (ndigits_fractional > 0) {
        val += fractional_part * pow10_negative<FType>(ndigits_fractional);
      }
    } else if (n_digits_whole == 0) {
      throw("neither whole part nor period");
    }
    // Consume exponent
    if (*p == 'e' || *p == 'E') {
      ++p;
      int ndigits_exponent;
      int exponent = parse_signed_integer(p, ndigits_exponent);
      if (ndigits_exponent == 0) {
        throw("found exponent char but no exponent");
      } else if (exponent > kMaxExponent || exponent < -kMaxExponent) {
        throw("exponent out of range");
      } else if (exponent > 0) {
        val *= pow10_positive<FType>(exponent);
      } else if (exponent < 0) {
        val *= pow10_negative<FType>(-exponent);
      }
    }
    return val * sign;
  }
};

#define TEXT_PARSER_TYPE(_type, _convert, ...)                                 \
  template <> inline text_parser& text_parser::operator>><_type>(_type& val) { \
    skip_space();                                                              \
    val = _convert(str, ##__VA_ARGS__);                                        \
    return *this;                                                              \
  }

TEXT_PARSER_TYPE(float, parse_fp<float>)
TEXT_PARSER_TYPE(double, parse_fp<double>)
TEXT_PARSER_TYPE(short, parse_signed_integer)
TEXT_PARSER_TYPE(int, parse_signed_integer)
TEXT_PARSER_TYPE(long, parse_signed_integer)
TEXT_PARSER_TYPE(long long, parse_signed_integer)
TEXT_PARSER_TYPE(unsigned short, parse_integer)
TEXT_PARSER_TYPE(unsigned int, parse_integer)
TEXT_PARSER_TYPE(unsigned long, parse_integer)
TEXT_PARSER_TYPE(unsigned long long, parse_integer)

#undef TEXT_PARSER_TYPE

}  // util
