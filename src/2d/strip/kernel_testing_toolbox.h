#ifndef KERNEL_TESTING_TOOLBOX
#define KERNEL_TESTING_TOOLBOX

#include <cmath>

namespace PET2D {

namespace Strip {
namespace Testing {

template <typename F> class FrameEvent;

template <typename F> struct Event {

  Event(F x, F y, F theta) : tan(std::tan(theta)), theta(theta), x(x), y(y){};
  Event(FrameEvent<F> fe, F R)
      : tan((fe.zup - fe.zdn) / (F(2.0) * R)),
        theta(atan(tan)),
        y(-F(0.5) * fe.dl * std::cos(theta)),
        x(F(0.5) * (fe.zup + fe.zdn + 2 * y * tan)){};

  const F tan;
  const F theta;
  const F y;
  const F x;
};

template <typename F> struct FrameEvent {

  FrameEvent(F zup, F zdn, F dl) : zup(zup), zdn(zdn), dl(dl) {}
  FrameEvent(Event<F> evt, F R)
      : zup(evt.x + (R - evt.y) * evt.tan),
        zdn(evt.x - (R + evt.y) * evt.tan),
        dl(-F(2) * evt.y * std::sqrt(F(1) + evt.tan * evt.tan)){};

  const F zup, zdn, dl;
};
}
}
}

#endif  // KERNEL_TESTING_TOOLBOX
