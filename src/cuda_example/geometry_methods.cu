#ifndef _GEOMETRY_METHODS_H_
#define _GEOMETRY_METHODS_H_

#include "data_structures.cu"
#include "config.h"
//#include "point.cu"

CUDA_CALLABLE_MEMBER Secant_Points
secant(float x, float y, float angle, float radius) {

  float a = std::sin(angle);
  float b = -std::cos(angle);
  float c = a * x + b * y;

  // std::cout << "a: " << a << " b: " << b << " c: " << c << std::endl;

  // helper variables
  float b2 = b * b;
  float b2c = b2 * c;
  float ac = a * c;
  float a2_b2 = a * a + b2;
  float b_a2_b2 = b * a2_b2;

  float sq = sqrt(b2 * (-(c * c) + a2_b2 * radius * radius));
  float asq = a * sq;

  // std::cout << "sq: " << sq << " asq: " << asq << std::endl;

  Secant_Points secant_positions;

  secant_positions.x1 = (ac - sq) / a2_b2;
  secant_positions.y1 = ((b2c + asq) / b_a2_b2);
  secant_positions.x2 = (ac + sq) / a2_b2;
  secant_positions.y2 = ((b2c - asq) / b_a2_b2);

  return secant_positions;
}

CUDA_CALLABLE_MEMBER Secant_Angle secant_angles(Secant_Points& e) {

  Secant_Angle temp;
  temp.angle1 = atan2(e.y1, e.x1);
  temp.angle2 = atan2(e.y2, e.x2);

  return temp;
}

CUDA_CALLABLE_MEMBER int section(float angle, int n_detectors) {
  // converting angles to [0,2 Pi) interval
  float normalised_angle = angle > 0 ? angle : (float)2.0 * M_PI + angle;
  return static_cast<int>(round(normalised_angle * n_detectors * INV_TWO_PI)) %
         (n_detectors);
}

CUDA_CALLABLE_MEMBER Secant_Sections
secant_sections(Secant_Points& e, int n_detectors) {

  Secant_Angle angles = secant_angles(e);

  Secant_Sections temp;

  temp.ss1 = section(angles.angle1, n_detectors);
  temp.ss2 = section(angles.angle2, n_detectors);

  return temp;
}

CUDA_CALLABLE_MEMBER int intersections(float x,
                                       float y,
                                       float angle,
                                       Detector_Ring& ring,
                                       int detector_id,
                                       Hits& hit) {

  // float temp_x =
  // ring.detector_list[detector_id].get_element(3).get_locationx();
  // float temp_y =
  // ring.detector_list[detector_id].get_element(3).get_locationy();

  float p1_x = ring.detector_list[detector_id].points[3].x;
  float p1_y = ring.detector_list[detector_id].points[3].y;

  // std::cout << "temp_x: " <<  temp_x << " " << "temp_y: " << temp_y <<
  // std::endl;

  float a = std::sin(angle);
  float b = -std::cos(angle);
  float c = a * x + b * y;

  // std::cout << "a: " << a << " " << "b: " << b << " " << "c: " << c <<
  // std::endl;

  float v1 = a * p1_x + b * p1_y - c;

  // std::cout <<"V1: " <<  v1 << std::endl;

  int r = 0;

  for (int i = 0; i < 4; i++) {

    float p2_x = ring.detector_list[detector_id].points[i].x;
    float p2_y = ring.detector_list[detector_id].points[i].y;

    float v2 = a * p2_x + b * p2_y - c;

    if (v2 == 0.0f) {
      hit.p[r].x = ring.detector_list[detector_id].points[i].x;
      hit.p[r].y = ring.detector_list[detector_id].points[i].y;
      // v2 is crossing point
      r++;
      // std::cout << " " << p2_x << "  " << p2_y;
      if (r == 2)
        return r;
    } else if (v1 * v2 < 0.0f) {
      // calculate intersection

      float m = a * (p1_x - p2_x) + b * (p1_y - p2_y);
      /*
       std::cout
       << (c * (p1_x - p2_x) + b * (p2_x * p1_y - p1_x * p2_y)) / m
       << "  "
       << ((c * (p1_y - p2_y) + a * (p1_x * p2_y - p2_x * p1_y))
       / m) << std::endl;*/
      hit.p[r].x = (c * (p1_x - p2_x) + b * (p2_x * p1_y - p1_x * p2_y)) / m;
      hit.p[r].y = (c * (p1_y - p2_y) + a * (p1_x * p2_y - p2_x * p1_y)) / m;

      r++;

      if (r == 2)
        return r;
    }
    v1 = v2;
    p1_x = p2_x;
    p1_y = p2_y;
  }
  return r;
}

CUDA_CALLABLE_MEMBER float secant_angle(float x1, float y1) {
  return atan2(y1, x1);
}

CUDA_CALLABLE_MEMBER bool check_for_hits(int inner,
                                         int outer,
                                         int x,
                                         int y,
                                         float angle,
                                         int n_detectors,
                                         Detector_Ring& ring,
                                         int& detector,
                                         Hits& hit) {

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

#endif
