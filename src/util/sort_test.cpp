#include "test.h"
#include "sort.h"
#include "array.h"

TEST("util/sort") {

  SECTION("simple_sort") {
    util::array<4, int> to_sort(3, 2, 4, 1);
    util::array<4, int> sorted(1, 2, 3, 4);
    util::heap_sort(
        to_sort.begin(), to_sort.end(), [](int a, int b) { return a < b; });
    REQUIRE(to_sort == sorted);
  }

  SECTION("index_sort") {
    util::array<4, double> distance(2., .2, 4., .1);
    util::array<4, int> to_sort(0, 1, 2, 3);
    util::array<4, int> sorted(3, 1, 0, 2);
    util::heap_sort(to_sort.begin(),
                    to_sort.end(),
                    [&](int a, int b) { return distance[a] < distance[b]; });
    REQUIRE(to_sort == sorted);
  }
}
