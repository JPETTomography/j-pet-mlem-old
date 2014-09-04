#ifndef THREAD_SAFE_STATS_H
#define THREAD_SAFE_STATS_H

#define DEFINE_STAT(T, name) \
  std::vector<T> name;       \
  T total_##name;

#define DEFINE_STAT_INITIALIZER(name) name(n_threads_, 0), total_##name(0)

#define FILL_STAT_WITH(name, value) std::fill_n(name.begin(), n_threads_, value)
#define COLLECT_STAT(name) total_##name += name[i];

template <typename T> struct Stats {
  Stats(int n_threads)
      : n_threads_(n_threads),
        DEFINE_STAT_INITIALIZER(n_events_processed_),
        DEFINE_STAT_INITIALIZER(n_pixels_processed_),
        DEFINE_STAT_INITIALIZER(n_kernel_calls_),
        DEFINE_STAT_INITIALIZER(bb_width_sum_),
        DEFINE_STAT_INITIALIZER(bb_height_sum_),
        DEFINE_STAT_INITIALIZER(bb_width2_sum_),
        DEFINE_STAT_INITIALIZER(bb_height2_sum_),
        DEFINE_STAT_INITIALIZER(bb_width_height_sum_) {}

  int n_threads_;
  DEFINE_STAT(T, n_events_processed_)
  DEFINE_STAT(T, n_pixels_processed_)
  DEFINE_STAT(T, n_kernel_calls_)
  DEFINE_STAT(T, bb_width_sum_)
  DEFINE_STAT(T, bb_height_sum_)
  DEFINE_STAT(T, bb_width2_sum_)
  DEFINE_STAT(T, bb_height2_sum_)
  DEFINE_STAT(T, bb_width_height_sum_)

  void fill(T value) {
    FILL_STAT_WITH(n_events_processed_, value);
    FILL_STAT_WITH(n_pixels_processed_, value);
    FILL_STAT_WITH(n_kernel_calls_, value);
    FILL_STAT_WITH(bb_width_sum_, value);
    FILL_STAT_WITH(bb_height_sum_, value);
    FILL_STAT_WITH(bb_width2_sum_, value);
    FILL_STAT_WITH(bb_height2_sum_, value);
    FILL_STAT_WITH(bb_width_height_sum_, value);
  }

  void collect() {
    for (int i = 0; i < n_threads_; i++) {
      COLLECT_STAT(n_events_processed_);
      COLLECT_STAT(n_pixels_processed_);
      COLLECT_STAT(n_kernel_calls_);
      COLLECT_STAT(bb_width_sum_)
      COLLECT_STAT(bb_height_sum_)
      COLLECT_STAT(bb_width2_sum_)
      COLLECT_STAT(bb_height2_sum_)
      COLLECT_STAT(bb_width_height_sum_) {}
    }
  }
};

#endif  // THREAD_SAFE_STATS_H
