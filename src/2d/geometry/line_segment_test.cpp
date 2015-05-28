#include "util/test.h"

#include "line_segment.h"

TEST("PET2D/geometry/line_segment") {
    using Point = PET2D::Point<float>;

    Point start(0,0);
    Point end(1,2);

    PET2D::LineSegment<float> segment(start, end);


}
