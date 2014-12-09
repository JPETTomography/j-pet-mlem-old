#include "cuda/compat.h"

namespace util {

namespace {
template <typename T, class Comparator>
_ void heap_sort_sift_down(T* a, int root, int count, Comparator compare) {
  // traverse down
  for (int child = 2 * root + 1; child < count; child = 2 * root + 1) {
    // select largest child of two
    if (child + 1 < count && compare(a[child], a[child + 1])) {
      ++child;
    }
    // if child is larger than parent, swap them
    if (compare(a[root], a[child])) {
      compat::swap(a[child], a[root]);
      root = child;
    } else {
      break;
    }
  }
}
}

/// Sorts given data using \c begin and \c end iterators

/// This implements heapsort was invented by J. W. J. Williams with
/// \f$ O(n \log n) \f$ worst and \f$ \Omega(n), O(n \log n) \f$ best case
/// performance.
///
/// This is somehow drop-in repleacement for \c std::sort that is compatible
/// with CUDA and requires no additional memory (sorts in-place).
template <typename T, class Comparator>
_ void heap_sort(T* begin,           ///< data begin iterator
                 T* end,             ///< data end iterator
                 Comparator compare  ///< comparator
                 ) {
  int count = static_cast<int>(end - begin);
  // heapify
  for (int root = (count - 2) / 2; root >= 0; --root) {
    heap_sort_sift_down(begin, root, count, compare);
  }
  // sort
  for (int last = count - 1; last > 0; --last) {
    compat::swap(begin[last], begin[0]);  // move largest value to the end
    heap_sort_sift_down(begin, 0, last, compare);
  }
}
}  // util
