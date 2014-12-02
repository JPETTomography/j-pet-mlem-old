#include "cuda/compat.h"

namespace util {

namespace {
template <typename T, class Comparator>
_ void heap_sort_sift_down(T* a, int start, int end, Comparator compare) {
  int root = start;
  while (root * 2 + 1 < end) {
    int child = 2 * root + 1;
    if (child + 1 < end && compare(a[child], a[child + 1])) {
      ++child;
    }
    if (compare(a[root], a[child])) {
      compat::swap(a[child], a[root]);
      root = child;
    } else {
      return;
    }
  }
}
}

template <typename T, class Comparator>
_ void heap_sort(T* a, T* b, Comparator compare) {
  int count = static_cast<int>(b - a);
  // heapify
  for (int start = (count - 2) / 2; start >= 0; --start) {
    heap_sort_sift_down(a, start, count, compare);
  }
  // sort
  for (int end = count - 1; end > 0; --end) {
    compat::swap(a[end], a[0]);
    heap_sort_sift_down(a, 0, end, compare);
  }
}
}  // util
