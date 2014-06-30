#pragma once

#include <cstdio>
#include <cmath>

// detect build variant
#if _OPENMP && HAVE_CUDA
#define VARIANT "OpenMP/CUDA"
#elif _OPENMP
#define VARIANT "OpenMP"
#elif HAVE_CUDA
#define VARIANT "CUDA"
#else
#define VARIANT "single-threaded CPU"
#endif

#if _OPENMP
#include <omp.h>
#endif

class Progress {
 public:
  Progress(bool enable,
           unsigned long total,
           unsigned long reasonable_update =
               std::numeric_limits<unsigned long>::max())
      : enable(enable),
        total(total),
        start_time(time(NULL)),
        mask(1),
        last_completed(std::numeric_limits<unsigned long>::max()) {
    // computes mask that shows percentage only ~ once per thousand of total
    auto total_resonable_update = total / 1000;
    while (mask < total_resonable_update && mask < reasonable_update) {
      mask <<= 1;
    }
    --mask;
  }

  void operator()(unsigned long completed) {

    // limit updates so they are not too often
    if (!enable || (completed & mask) != 0 || last_completed == completed)
      return;
#if _OPENMP
    if (omp_get_thread_num() != 0)
      return;
#endif

    last_completed = completed;

    double persec = (double)completed / (double)(time(NULL) - start_time);

    std::cerr << " " << std::round((double)completed / (double)total * 100.0)
              << "% " << completed << "/" << total;

    if (!std::isnan(persec)) {
      std::cerr << " "
                << timetostr(std::round((double)(total - completed) / persec))
                << " left, ";
      std::cerr << timetostr(std::round((double)total / persec)) << " total";
    }
    std::cerr << "\r";
  }

 private:
  static const char* timetostr(int sec) {
    static char out[64];
    int min = sec / 60;
    sec = sec % 60;
    int hour = min / 60;
    min = min % 60;
    std::sprintf(out, "%2d:%02d:%02d", hour, min, sec);
    return out;
  }

  bool enable;
  unsigned long total;
  time_t start_time;
  unsigned long mask;
  unsigned long last_completed;
};
