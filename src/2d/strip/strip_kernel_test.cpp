#include "util/test.h"

#include "common/types.h"

#include "2d/strip/kernel_testing_toolbox.h"

using namespace PET2D::Strip::Testing;

TEST("strip/events") {

  Event<F> evt(0.1, -.2, M_PI / 3);

  CHECK(evt.x == Approx(0.1));
  CHECK(evt.y == Approx(-0.2));
  CHECK(evt.theta == Approx(M_PI / 3));
  CHECK(evt.tan == Approx(tan(M_PI / 3)));

  FrameEvent<F> fevnt(-.2, .1, 0.3);

  CHECK(fevnt.zup == Approx(-0.2));
  CHECK(fevnt.zdn == Approx(0.1));
  CHECK(fevnt.dl == Approx(0.3));
}

TEST("strip/event/conversion/Event-FrameEvent-Event") {
  F R = 0.4;
  Event<F> evt(0.1, -0.2, M_PI / 6);

  Event<F> conv(FrameEvent<F>(evt, R), R);

  CHECK(conv.x == Approx(evt.x));
  CHECK(conv.y == Approx(evt.y));
  CHECK(conv.theta == Approx(evt.theta));
  CHECK(conv.tan == Approx(evt.tan));
}
