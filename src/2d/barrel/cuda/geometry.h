#pragma once

#include "config.h"

namespace PET2D {
namespace Barrel {

/// GPU internal namespace

/// \todo This namespace should disappear, and we shall use shared classes like
/// it is done for PET2D::Strip.
namespace GPU {

/// \cond PRIVATE

struct Point {
  float x, y;
};

struct Hits {
  Point p[2];
};

struct Detector {
  Point points[4];
};

struct DetectorRing {
  Detector detector_list[NUMBER_OF_DETECTORS];
};

struct MatrixElement {
  float hit[LORS];
};

struct SecantPoints {
  float x1, y1, x2, y2;
};

struct SecantAngles {
  float angle1, angle2;
};

struct SecantSections {
  int ss1, ss2;
};

struct LOR {
  int lor_a;
  int lor_b;
  int index() const { return (lor_a * (lor_a + 1)) / 2 + lor_b; }
};

/// \endcond

}  // GPU
}  // Barrel
}  // PET2D
