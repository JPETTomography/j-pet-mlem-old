#pragma once

#include "config.h"

struct data {

  float x, y;
  float f_data[8];
};

// x,y,z,w
//__builtin_align__(32)

struct Points {

  float x, y;
};

struct Hits {

  Points p[2];
};

struct Detectors {

  Points points[4];
};

struct Detector_Ring {

  Detectors detector_list[NUMBER_OF_DETECTORS];
};

struct Matrix_Element {

  float lor[LORS];
};

struct Secant_Points {

  float x1, y1, x2, y2;
};

struct Secant_Angle {

  float angle1, angle2;
};

struct Secant_Sections {

  int ss1, ss2;
};
