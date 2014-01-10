#pragma once

#include <cstdio>

static inline const char* report_time(int sec) {
  static char out[64];
  int min = sec / 60;
  sec = sec % 60;
  int hour = min / 60;
  min = min % 60;
  std::sprintf(out, "%2d:%02d:%02d", hour, min, sec);
  return out;
}

inline void report_progress(time_t start_time, int completed, int total) {

  double persec = (double)completed / (double)(time(NULL) - start_time);

  std::cerr << " " << std::round((double)completed / (double)total * 100.0)
            << "% " << completed << "/" << total;

  if (!std::isnan(persec)) {
    std::cerr << " "
              << report_time(std::round((double)(total - completed) / persec))
              << " left, ";
    std::cerr << report_time(std::round((double)total / persec)) << " total";
  }
  std::cerr << "\r";
}
