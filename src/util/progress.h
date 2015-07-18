#pragma once

#include <iostream>
#include <ostream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <chrono>

#if _OPENMP
#include <omp.h>
#endif

namespace util {

/// Provides automatic progress and \a EST display.
class progress {
 public:
  typedef std::chrono::high_resolution_clock clock_t;
  typedef clock_t::time_point time_t;

  enum {
    Disabled = 0,  ///< Do not display any progress.
    Estimate,      ///< Show only one line progress with estimated time.
    Benchmark,     ///< Show detailed statistics per each update.
  };

  /// Constructs new progress handler.
  progress(
      int verbosity,                     ///< Verbosity level (0-2)
      unsigned long total,               ///< Total number of iterations
      unsigned long reasonable_update =  ///< Iterations that trigger update
      std::numeric_limits<unsigned long>::max())
      : verbosity(verbosity),
        total(total),
        total_m_1(total - 1),
        start_time(clock_t::now()),
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

  /// Constructs new progress handler.
  progress(
      bool enabled,                      ///< If progress is enabled at all
      unsigned long total,               ///< Total number of iterations
      unsigned long reasonable_update =  ///< Iterations that trigger update
      std::numeric_limits<unsigned long>::max())
      : progress(enabled ? Estimate : Disabled, total, reasonable_update) {}

  /// Call when begining or finishing iteration
  void operator()(unsigned long completed,  ///< Number of completed so far
                  bool finished = false     ///< Whether we finish iteration
                  ) {

    // limit updates so they are not too often
    if (!(verbosity &&
          ((completed == total_m_1 && finished) || (completed & mask) == 0)))
      return;
#if _OPENMP
    if (omp_get_thread_num() != 0)
      return;
#endif

    // Verbose (benchmark) mode
    if (verbosity >= Benchmark) {
      if (!finished) {
        start_time = clock_t::now();
      } else {
        auto ellapsed_ms = ellapsed() * 1000;
        total_ellapsed_ms += ellapsed_ms;
        auto prev_precision = std::cout.precision();
        std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(2)
                  << std::setw(4) << (completed + 1) << " " << std::setw(10)
                  << ellapsed_ms << std::endl;
        if (completed == total_m_1) {
          std::cout << "# av " << std::setw(10) << total_ellapsed_ms / total
                    << "   total " << total_ellapsed_ms << std::endl;
        }
        std::cout << std::resetiosflags(std::ios::fixed)
                  << std::setprecision(prev_precision);
      }
      return;
    }
    // Estimate time mode
    else {
      if (finished) {
        if (completed < total_m_1)
          return;
        ++completed;
      }

      if (completed == last_completed)
        return;
      last_completed = completed;

      double persec = (double)completed / ellapsed();

      std::cerr << " " << std::round((double)completed / total * 100.0) << "% "
                << completed << "/" << total;

      if (!std::isnan(persec) && completed > 0) {
        std::cerr << " "
                  << timetostr(std::round((double)(total - completed) / persec))
                  << " left, ";
        std::cerr << timetostr(std::round((double)total / persec))
                  << " total  ";
      }
      if (finished && completed == total) {
        std::cerr << std::endl;
      } else {
        std::cerr << "\r";
      }
    }
  }

 private:
  double ellapsed() {
    return std::chrono::duration_cast<std::chrono::duration<double>>(
               clock_t::now() - start_time)
        .count();
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
  unsigned long total_m_1;
  time_t start_time;
  unsigned long mask;
  unsigned long last_completed;
  double total_ellapsed_ms;
};
}  // util
