#include <ostream>

// this backtrace information works only on Linux && OS X
#if (__linux || __APPLE__) && (__clang__ || __GNUC__)

#include <string>
#include <memory>

#include <dlfcn.h>
#include <execinfo.h>
#include <typeinfo>
#include <cxxabi.h>

namespace {
void* last_frames[20];
size_t last_size;
}

extern "C" {
void __cxa_throw(void* ex,
#if __clang__
                 std::type_info* info,
#else
                 void* info,
#endif
                 void (*dest)(void*)) {
  last_size =
      backtrace(last_frames, sizeof(last_frames) / sizeof(*last_frames));
  static void (*const rethrow)(void*, void*, void (*)(void*))
      __attribute__((noreturn)) =
#if __clang__
          (__attribute__((noreturn)) void (*)(
              void*, void*, void (*)(void*)))dlsym(RTLD_NEXT, "__cxa_throw");
#else
          (void (*)(void*, void*, void (*)(void*)))dlsym(RTLD_NEXT,
                                                         "__cxa_throw");
#endif
  rethrow(ex, info, dest);
}
}

namespace util {

void print_backtrace(std::ostream& out) {
  char** symbols = backtrace_symbols(last_frames, last_size);
  for (size_t frame = 1; frame < last_size; ++frame) {
    std::string frame_str(symbols[frame]);
    Dl_info info;
    if (dladdr(last_frames[frame], &info) && info.dli_sname) {
      if (info.dli_sname[0] == '_') {
        int status = -1;
        const char* demangled =
            abi::__cxa_demangle(info.dli_sname, 0, 0, &status);
        if (demangled && !status) {
          std::string sname_str(info.dli_sname);
          size_t start_pos = frame_str.find(sname_str);
          if (start_pos != std::string::npos) {
            frame_str.replace(start_pos, sname_str.length(), demangled);
          }
        }
        if (demangled) {
          std::free((void*)demangled);
        }
      }
      out << frame_str << std::endl;
    }
  }
}

}  // util

#else
namespace util {
void print_backtrace(std::ostream& out) { (void)(out); }
}
#endif
