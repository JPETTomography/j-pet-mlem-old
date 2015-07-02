#include <ostream>

// http://stackoverflow.com/questions/19200207/redirect-debug-output-to-null

namespace util {

template <class cT, class traits = std::char_traits<cT>>
class basic_null_buf : public std::basic_streambuf<cT, traits> {
  typename traits::int_type overflow(typename traits::int_type c) {
    return traits::not_eof(c);  // indicate success
  }
};

template <class cT, class traits = std::char_traits<cT>>
class basic_null_ostream : public std::basic_ostream<cT, traits> {
 public:
  basic_null_ostream()
      : std::basic_ios<cT, traits>(&sbuf),
        std::basic_ostream<cT, traits>(&sbuf) {
    this->init(&sbuf);
  }

 private:
  basic_null_buf<cT, traits> sbuf;
};

using null_ostream = basic_null_ostream<char>;
using null_wostream = basic_null_ostream<wchar_t>;
}
