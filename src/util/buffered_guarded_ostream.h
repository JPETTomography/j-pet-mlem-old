#ifndef BUFFERED_GUARDED_OSTREAM
#define BUFFERED_GUARDED_OSTREAM

#include <mutex>
#include <vector>
#include <sstream>

template <typename S> class BufferedGuardedStream {
 public:
  BufferedGuardedStream(S& stream, int n_threads, int buffer_length)
      : n_threads(n_threads),
        stream_(stream),
        buffer_length_(buffer_length),
        buffer_(n_threads),
        index(n_threads) {
    for (int i = 0; i < n_threads; ++i) {
      buffer_[i].assign(buffer_length_, "");
      index[i] = 0;
    }
  };

  BufferedGuardedStream(S& stream, int n_threads)
      : BufferedGuardedStream(stream, n_threads, 512){};

  template <typename T> void writeln(int thread, const T& elem) {
    std::ostringstream ss;
    ss << elem << "\n";
    buffer_[thread][index[thread]] = ss.str();
    index[thread]++;
    if (index[thread] == buffer_length_) {
      {
        std::lock_guard<std::mutex> event_lock(event_mutex);
        for (int i = 0; i < index[thread]; ++i)
          stream_ << buffer_[thread][i];
        index[thread] = 0;
      }
    }
  }

  void flush() {
    {

      for (int thread = 0; thread < n_threads; ++thread) {
        for (int i = 0; i < index[thread]; ++i) {
          stream_ << buffer_[thread][i];
        }
        index[thread] = 0;
      }
    }
  }

  const int n_threads;

 private:
  std::mutex event_mutex;
  S& stream_;
  int buffer_length_;
  std::vector<std::vector<std::string>> buffer_;
  std::vector<int> index;
};

#endif  // BUFFERED_GUARDED_OSTREAM
