#pragma once

namespace PET2D {
namespace Strip {

/// Gathers various reconstruction statistics
template <typename T> struct ReconstructionStats {

#if USE_STATISTICS

#define DEFINE_STAT_INITIALIZER(name) thread_##name(n_threads, 0), name(0)

  Stats(size_t n_threads)
      : n_threads(n_threads),
        DEFINE_STAT_INITIALIZER(n_events_processed),
        DEFINE_STAT_INITIALIZER(n_pixels_processed),
        DEFINE_STAT_INITIALIZER(n_kernel_calls),
        DEFINE_STAT_INITIALIZER(bb_width_sum),
        DEFINE_STAT_INITIALIZER(bb_height_sum),
        DEFINE_STAT_INITIALIZER(bb_width2_sum),
        DEFINE_STAT_INITIALIZER(bb_height2_sum),
        DEFINE_STAT_INITIALIZER(bb_width_height_sum) {}

  size_t n_threads;

#define DEFINE_STAT(name)                                                    \
  std::vector<T> thread_##name;                                              \
  T name;                                                                    \
  void name##_by(size_t thread, T value) { thread_##name[thread] += value; } \
  void name##_by(T value = T(1)) { name##_by(omp_get_thread_num(), value); }

#else

  ReconstructionStats(size_t) {}

#define DEFINE_STAT(name)      \
  T name;                      \
  void name##_by(size_t, T) {} \
  void name##_by(T = T(1)) {}

#endif

  DEFINE_STAT(n_events_processed)
  DEFINE_STAT(n_pixels_processed)
  DEFINE_STAT(n_kernel_calls)
  DEFINE_STAT(bb_width_sum)
  DEFINE_STAT(bb_height_sum)
  DEFINE_STAT(bb_width2_sum)
  DEFINE_STAT(bb_height2_sum)
  DEFINE_STAT(bb_width_height_sum)

#if USE_STATISTICS

#define FILL_STAT_WITH(name, value) \
  std::fill_n(thread_##name.begin(), n_threads, value)

  void fill(T value = T()) {
    FILL_STAT_WITH(n_events_processed, value);
    FILL_STAT_WITH(n_pixels_processed, value);
    FILL_STAT_WITH(n_kernel_calls, value);
    FILL_STAT_WITH(bb_width_sum, value);
    FILL_STAT_WITH(bb_height_sum, value);
    FILL_STAT_WITH(bb_width2_sum, value);
    FILL_STAT_WITH(bb_height2_sum, value);
    FILL_STAT_WITH(bb_width_height_sum, value);
  }

#define COLLECT_STAT(name) name += thread_##name[i]

  void collect() {
    for (size_t i = 0; i < n_threads; i++) {
      COLLECT_STAT(n_events_processed);
      COLLECT_STAT(n_pixels_processed);
      COLLECT_STAT(n_kernel_calls);
      COLLECT_STAT(bb_width_sum);
      COLLECT_STAT(bb_height_sum);
      COLLECT_STAT(bb_width2_sum);
      COLLECT_STAT(bb_height2_sum);
      COLLECT_STAT(bb_width_height_sum);
    }
  }

#else

  void fill(T = T()) {}
  void collect() {}

#endif
};
}  // Strip
}  // PET2D
