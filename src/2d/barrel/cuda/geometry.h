#pragma once

#include "config.h"
#include "2d/geometry/point.h"
#include "2d/barrel/lor.h"
#include "util/cuda/compat.h"
#include "prng.cuh"

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

_ inline void create_detector_ring(int i,
                                   float h_detector,
                                   float w_detector,
                                   float radius,
                                   DetectorRing& dr) {

  SquareDetector detector_base;

  detector_base.points[0].x = (h_detector / 2) + radius + (h_detector / 2);
  detector_base.points[0].y = w_detector / 2;
  detector_base.points[1].x = (h_detector / 2) + radius + (h_detector / 2);
  detector_base.points[1].y = -w_detector / 2;
  detector_base.points[2].x = (-h_detector / 2) + radius + (h_detector / 2);
  detector_base.points[2].y = -w_detector / 2;
  detector_base.points[3].x = (-h_detector / 2) + radius + (h_detector / 2);
  detector_base.points[3].y = w_detector / 2;

  dr.detector_list[i] = detector_base;

  float angle = 2 * (float)M_PI * i / NUMBER_OF_DETECTORS;

  float sin_phi = compat::sin(angle);
  float cos_phi = compat::cos(angle);

  for (int j = 0; j < 4; ++j) {

    float temp_x = dr.detector_list[i].points[j].x;
    float temp_y = dr.detector_list[i].points[j].y;

    dr.detector_list[i].points[j].x = temp_x * cos_phi - temp_y * sin_phi;
    dr.detector_list[i].points[j].y = temp_x * sin_phi + temp_y * cos_phi;
  }
}

_ inline void secant(SecantPoints& s1,
                     SecantPoints& s2,
                     float x,
                     float y,
                     float angle,
                     float radius_s1,
                     float radius_s2) {

  float a = compat::sin(angle);
  float b = -compat::cos(angle);
  float c = a * x + b * y;

  float b2 = b * b;
  float b2c = b2 * c;
  float ac = a * c;
  float a2_b2 = a * a + b2;
  float b_a2_b2 = b * a2_b2;

  float sq_s1 = sqrt(b2 * (-(c * c) + a2_b2 * radius_s1 * radius_s1));
  float sq_s2 = sqrt(b2 * (-(c * c) + a2_b2 * radius_s2 * radius_s2));

  float asq_s1 = a * sq_s1;
  float asq_s2 = a * sq_s2;

  s1.x1 = (ac - sq_s1) / a2_b2;
  s1.y1 = ((b2c + asq_s1) / b_a2_b2);
  s1.x2 = (ac + sq_s1) / a2_b2;
  s1.y2 = ((b2c - asq_s1) / b_a2_b2);

  s2.x1 = (ac - sq_s2) / a2_b2;
  s2.y1 = ((b2c + asq_s2) / b_a2_b2);
  s2.x2 = (ac + sq_s2) / a2_b2;
  s2.y2 = ((b2c - asq_s2) / b_a2_b2);
}

_ inline SecantAngles secant_angles(SecantPoints& e) {

  SecantAngles temp;
  temp.angle1 = atan2(e.y1, e.x1);
  temp.angle2 = atan2(e.y2, e.x2);

  return temp;
}

_ inline int section(float angle, int n_detectors) {
  // converting angles to [0,2 Pi) interval
  float normalised_angle = angle > 0 ? angle : 2 * (float)M_PI + angle;
  return static_cast<int>(round(normalised_angle * n_detectors * 0.1591549f)) %
         (n_detectors);
}

_ inline SecantSections secant_sections(SecantPoints& e, int n_detectors) {

  SecantAngles angles = secant_angles(e);

  SecantSections temp;

  temp.ss1 = section(angles.angle1, n_detectors);
  temp.ss2 = section(angles.angle2, n_detectors);

  return temp;
}

_ inline int intersections(float x,
                           float y,
                           float angle,
                           DetectorRing& ring,
                           int& detector_id,
                           Hits& hit) {

  float p1_x = ring.detector_list[detector_id].points[3].x;
  float p1_y = ring.detector_list[detector_id].points[3].y;

  float a = compat::sin(angle);
  float b = -compat::cos(angle);
  float c = a * x + b * y;

  float v1 = a * p1_x + b * p1_y - c;

  int r = 0;

  for (int i = 0; i < 4; i++) {

    float p2_x = ring.detector_list[detector_id].points[i].x;
    float p2_y = ring.detector_list[detector_id].points[i].y;

    float v2 = a * p2_x + b * p2_y - c;

    if (v2 == 0) {
      hit.p[r].x = ring.detector_list[detector_id].points[i].x;
      hit.p[r].y = ring.detector_list[detector_id].points[i].y;

      r++;

      if (r == 2) {
        return r;
      }
    } else if (v1 * v2 < 0) {
      // calculate intersection

      float m = a * (p1_x - p2_x) + b * (p1_y - p2_y);

      hit.p[r].x = (c * (p1_x - p2_x) + b * (p2_x * p1_y - p1_x * p2_y)) / m;
      hit.p[r].y = (c * (p1_y - p2_y) + a * (p1_x * p2_y - p2_x * p1_y)) / m;

      r++;

      if (r == 2) {
        return r;
      }
    }
    v1 = v2;
    p1_x = p2_x;
    p1_y = p2_y;
  }
  return r;
}

_ inline float SecantAngles(float x1, float y1) { return atan2(y1, x1); }

_ inline float length(Point p1, Point p2) {

  return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}

_ inline float nearest_distance(Point p1, Point p2, Point center) {

  return compat::min(length(p1, center), length(p2, center));
}

_ inline bool check_for_hits(int flag,
                             int inner,
                             int outer,
                             float x,
                             float y,
                             float angle,
                             int n_detectors,
                             DetectorRing& ring,
                             int& detector,
                             Hits& hit,
                             unsigned int* seed,
                             float& depth) {

  if (!flag) {
    return false;
  }

  int points;

  int step = ((n_detectors + inner - outer) % (n_detectors) >
              (n_detectors + outer - inner) % (n_detectors))
                 ? 1
                 : n_detectors - 1;
  int end = (outer + step) % (n_detectors);
  for (int i = inner; i != end; i = (i + step) % (n_detectors)) {
    points = intersections(x, y, angle, ring, i, hit);

#if CLOCK_TEST
    exec_inter++;
#endif
    if (points == 2) {

      detector = i;

      depth =
          -compat::log(HybridTaus(seed[0], seed[1], seed[2], seed[3])) * 0.1f;

      if (depth < length(hit.p[0], hit.p[1])) {
        return true;
      }
    }
  }

  return false;
}

_ inline int quantize_position(float position,
                               float step_size,
                               float n_positions) {
  if (position < 0)
    return n_positions / 2 - 1.0f -
           static_cast<int>(floor(-position / step_size));
  else
    return static_cast<int>(floor(position / step_size)) + n_positions / 2;
}

_ inline float max_dl(float max_bias_size, float c_outer_radius) {
  return 2 * c_outer_radius + max_bias_size;
}

_ inline int n_positions(float step_size,
                         float max_bias_size,
                         float c_outer_radius) {
  // since position needs to be symmetric against (0,0) number must be even
  return ((int)(ceil(2 * max_dl(max_bias_size, c_outer_radius) / step_size)) +
          1) /
         2 * 2;
}

/// \endcond

}  // GPU
}  // Barrel
}  // PET2D
