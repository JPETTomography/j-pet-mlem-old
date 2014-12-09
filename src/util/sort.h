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

template <typename T, class Comparator>
_ void heap_sort(T* a, T* b, Comparator compare) {
  int count = static_cast<int>(b - a);
  // heapify
  for (int root = (count - 2) / 2; root >= 0; --root) {
    heap_sort_sift_down(a, root, count, compare);
  }
  // sort
  for (int end = count - 1; end > 0; --end) {
    compat::swap(a[end], a[0]);  // move largest value to the end
    heap_sort_sift_down(a, 0, end, compare);
  }
}
}  // util
