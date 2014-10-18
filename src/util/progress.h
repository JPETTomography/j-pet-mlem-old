#pragma once

#include <ostream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <chrono>

#if _OPENMP
#include <omp.h>
#endif

/// Provides automatic progress and \a EST display
class Progress {
 public:
  typedef std::chrono::high_resolution_clock Clock;
  typedef Clock::time_point Time;

  enum {
    Disabled = 0,
    Estimate,
    Benchmark,
  };

  Progress(int verbosity,
           unsigned long total,
           unsigned long reasonable_update =
               std::numeric_limits<unsigned long>::max())
      : verbosity(verbosity),
        total(total),
        start_time(Clock::now()),
        mask(1),
        last_completed(std::numeric_limits<unsigned long>::max()),
        total_ellapsed_ms(0) {
    // computes mask that shows percentage only ~ once per thousand of total
    auto total_resonable_update = total / 1000;
    while (mask < total_resonable_update && mask < reasonable_update) {
      mask <<= 1;
    }
    --mask;

    if (verbosity >= Benchmark) {
      std::cout << "# it     time (ms)" << std::endl;
    }
  }

  Progress(bool enabled,
           unsigned long total,
           unsigned long reasonable_update =
               std::numeric_limits<unsigned long>::max())
      : Progress(enabled ? Estimate : Disabled, total, reasonable_update) {}

  void operator()(unsigned long completed, bool finished = false) {

    // limit updates so they are not too often
    if (!verbosity || (completed & mask) != 0 || last_completed == completed)
      return;
#if _OPENMP
    if (omp_get_thread_num() != 0)
      return;
#endif

    // Verbose (benchmark) mode
    if (verbosity >= Benchmark) {
      if (!finished) {
        start_time = Clock::now();
      } else {
        auto ellapsed_ms = ellapsed() * 1000;
        total_ellapsed_ms += ellapsed_ms;
        auto prev_precision = std::cout.precision();
        std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(2)
                  << std::setw(4) << (completed + 1) << " " << std::setw(10)
                  << ellapsed_ms << std::endl;
        if (completed == total - 1) {
          std::cout << "# av " << std::setw(10) << total_ellapsed_ms / total
                    << "   total " << total_ellapsed_ms << std::endl;
        }
        std::cout << std::resetiosflags(std::ios::fixed)
                  << std::setprecision(prev_precision);
      }
      return;
    }

    // Estimate time mode
    last_completed = completed;

    double persec = (double)completed / ellapsed();

    std::cerr << " " << std::round((double)completed / total * 100.0) << "% "
              << completed << "/" << total;

    if (!std::isnan(persec) && completed > 0) {
      std::cerr << " "
                << timetostr(std::round((double)(total - completed) / persec))
                << " left, ";
      std::cerr << timetostr(std::round((double)total / persec)) << " total";
    }
    std::cerr << "\r";
  }

 private:
  double ellapsed() {
    return std::chrono::duration_cast<std::chrono::duration<double>>(
               Clock::now() - start_time).count();
  }

  static const char* timetostr(int sec) {
    static char out[64];
    int min = sec / 60;
    sec = sec % 60;
    int hour = min / 60;
    min = min % 60;
    std::sprintf(out, "%2d:%02d:%02d", hour, min, sec);
    return out;
  }

  int verbosity;
  unsigned long total;
  Time start_time;
  unsigned long mask;
  unsigned long last_completed;
  double total_ellapsed_ms;
};
