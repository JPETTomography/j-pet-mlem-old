#include "util/test.h"

#include "frame.h"

using Point = PET3D::Point<float>;

TEST("PET3D/geometry/frame") {

    Frame<float> frame(Point(-0.205, 0.0, 0.0), Point(0.205, 0.0, 0.0));

}
