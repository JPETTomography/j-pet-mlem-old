#include <random>

#include "test.h"
#include "random.h"

TEST("util/random") {
  util::random::tausworthe rng;
  util::random::uniform_real_distribution<double> one;
  util::random::uniform_real_distribution<double> d99to100(99., 100.);

  CHECK(one.scale<util::random::tausworthe>() == Approx(2.328306436538696e-10));

  for (auto i = 0; i < 100; ++i) {
    auto r = one(rng);
    CHECK(r >= 0.);
    CHECK(r <= 1.);
  }

  for (auto i = 0; i < 100; ++i) {
    auto r = d99to100(rng);
    CHECK(r >= 99.);
    CHECK(r <= 100.);
  }
}

TEST("util/random/restore") {
  std::random_device rd;
  util::random::tausworthe rng1(rd());
  util::random::tausworthe::state_type state1;
  rng1.save(state1);
  util::random::tausworthe rng2(state1);
  util::random::tausworthe::state_type state2;
  rng2.save(state2);
  util::random::tausworthe rng3(state2);
  CHECK(rng1 == rng2);
  CHECK(rng2 == rng3);
}
