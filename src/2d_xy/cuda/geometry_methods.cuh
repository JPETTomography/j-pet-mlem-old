#pragma once

#include <cuda_runtime.h>

#include "geometry.h"
#include "config.h"

using namespace gpu;

__device__ void create_detector_ring(float& h_detector,
                                     float& w_detector,
                                     float& radius,
                                     DetectorRing& test_ring) {

  Detector detector_base;

  detector_base.points[0].x =
      (h_detector / 2.0f) + radius + (h_detector / 2.0f);
  detector_base.points[0].y = w_detector / 2.0f;
  detector_base.points[1].x =
      (h_detector / 2.0f) + radius + (h_detector / 2.0f);
  detector_base.points[1].y = -w_detector / 2.0f;
  detector_base.points[2].x =
      (-h_detector / 2.0f) + radius + (h_detector / 2.0f);
  detector_base.points[2].y = -w_detector / 2.0f;
  detector_base.points[3].x =
      (-h_detector / 2.0) + radius + (h_detector / 2.0f);
  detector_base.points[3].y = w_detector / 2.0f;

  test_ring.detector_list[threadIdx.x] = detector_base;

  float angle = 2.0f * M_PI * threadIdx.x / NUMBER_OF_DETECTORS;
  float sin_phi = __sinf(angle);
  float cos_phi = __cosf(angle);

  for (int j = 0; j < 4; ++j) {

    float temp_x = test_ring.detector_list[threadIdx.x].points[j].x;
    float temp_y = test_ring.detector_list[threadIdx.x].points[j].y;

    test_ring.detector_list[threadIdx.x].points[j].x =
        temp_x * cos_phi - temp_y * sin_phi;
    test_ring.detector_list[threadIdx.x].points[j].y =
        temp_x * sin_phi + temp_y * cos_phi;
  }
}

__device__ SecantPoints secant(float x, float y, float angle, float radius) {

  float a = std::sin(angle);
  float b = -std::cos(angle);
  float c = a * x + b * y;

  float b2 = b * b;
  float b2c = b2 * c;
  float ac = a * c;
  float a2_b2 = a * a + b2;
  float b_a2_b2 = b * a2_b2;

  float sq = sqrt(b2 * (-(c * c) + a2_b2 * radius * radius));
  float asq = a * sq;

  SecantPoints secant_positions;

  secant_positions.x1 = (ac - sq) / a2_b2;
  secant_positions.y1 = ((b2c + asq) / b_a2_b2);
  secant_positions.x2 = (ac + sq) / a2_b2;
  secant_positions.y2 = ((b2c - asq) / b_a2_b2);

  return secant_positions;
}

__device__ SecantAngles secant_angles(SecantPoints& e) {

  SecantAngles temp;
  temp.angle1 = atan2(e.y1, e.x1);
  temp.angle2 = atan2(e.y2, e.x2);

  return temp;
}

__device__ int section(float angle, int n_detectors) {
  // converting angles to [0,2 Pi) interval
  float normalised_angle = angle > 0 ? angle : (float)2.0 * M_PI + angle;
  return static_cast<int>(round(normalised_angle * n_detectors * 0.1591549f)) %
         (n_detectors);
}

__device__ SecantSections secant_sections(SecantPoints& e, int n_detectors) {

  SecantAngles angles = secant_angles(e);

  SecantSections temp;

  temp.ss1 = section(angles.angle1, n_detectors);
  temp.ss2 = section(angles.angle2, n_detectors);

  return temp;
}

__device__ int intersections(float x,
                             float y,
                             float angle,
                             DetectorRing& ring,
                             int detector_id,
                             Hits& hit) {

  float p1_x = ring.detector_list[detector_id].points[3].x;
  float p1_y = ring.detector_list[detector_id].points[3].y;

  float a = std::sin(angle);
  float b = -std::cos(angle);
  float c = a * x + b * y;

  float v1 = a * p1_x + b * p1_y - c;

  int r = 0;

  for (int i = 0; i < 4; i++) {

    float p2_x = ring.detector_list[detector_id].points[i].x;
    float p2_y = ring.detector_list[detector_id].points[i].y;

    float v2 = a * p2_x + b * p2_y - c;

    if (v2 == 0.0f) {
      hit.p[r].x = ring.detector_list[detector_id].points[i].x;
      hit.p[r].y = ring.detector_list[detector_id].points[i].y;

      r++;

      if (r == 2) {
        return r;
      }
    } else if (v1 * v2 < 0.0f) {
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

__device__ float SecantAngles(float x1, float y1) { return atan2(y1, x1); }

__device__ bool check_for_hits(int inner,
                               int outer,
                               float x,
                               float y,
                               float angle,
                               int n_detectors,
                               DetectorRing& ring,
                               int& detector,
                               Hits& hit,unsigned int *seed) {

  int points;

  int step = ((n_detectors + inner - outer) % n_detectors >
              (n_detectors + outer - inner) % n_detectors)
                 ? 1
                 : n_detectors - 1;
  int end = (outer + step) % n_detectors;
  for (int i = inner; i != end; i = (i + step) % n_detectors) {
    points = intersections(x, y, angle, ring, i, hit);

    if (points == 2) {

      detector = i;


//      if (-log(HybridTaus(seed[0], seed[1], seed[2], seed[3])) * inv_unit_prob_ <
//          (sqrt((hit.p[1].x - hit.p[0].x) * (hit.p[1].x - hit.p[0].x) +
//                (hit.p[1].y - hit.p[0].y) * (hit.p[1].y - hit.p[0].y)))) {
//      }


      return true;
    }

    // check if we got 2 point intersection
    // then test the model against these points distance
    // if (points.size() == 2) {
    //   auto deposition_depth = model.deposition_depth(gen);
    //   if (deposition_depth < (points[1] - points[0]).length()) {
    //    detector = i;
    //    depth = deposition_depth;
    //    return true;
    //  }
  }

  return false;
}
