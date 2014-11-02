#include "util/test.h"

#include "random.h"

TEST_CASE("random") {
  util::random::tausworthe gen;
  util::random::uniform_real_distribution<> one;
  util::random::uniform_real_distribution<> d99to100(99., 100.);

  CHECK(one.scale<util::random::tausworthe>() == Approx(2.328306436538696e-10));

  for (auto i = 0; i < 100; ++i) {
    auto r = one(gen);
    CHECK(r >= 0.);
    CHECK(r <= 1.);
  }

  for (auto i = 0; i < 100; ++i) {
    auto r = d99to100(gen);
    CHECK(r >= 99.);
    CHECK(r <= 100.);
  }
}
