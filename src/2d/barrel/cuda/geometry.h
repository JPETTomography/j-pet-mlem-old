#pragma once

#include "config.h"
#include "2d/geometry/point.h"
#include "2d/barrel/lor.h"

namespace PET2D {
namespace Barrel {

/// GPU internal namespace

/// \todo TODO: This namespace should disappear, and we shall use shared classes
/// like it is done for PET2D::Strip.
namespace GPU {

/// \cond PRIVATE

using Point = PET2D::Point<float>;
using LOR = PET2D::Barrel::LOR<>;

struct Hits {
  Point p[2];
};

struct SquareDetector {
  Point points[4];
};

struct DetectorRing {
  SquareDetector detector_list[NUMBER_OF_DETECTORS];
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

/// \endcond

}  // GPU
}  // Barrel
}  // PET2D
