#define CATCH_CONFIG_MAIN

#include <iostream>
#include <fstream>

#include "catch.hpp"
#include "reconstruction.h"


//TEST FOR EVENT (300,0,45.0) R = 450.0


struct E{

  double y;
  double z;
  double angle;

};

TEST_CASE("ovec", "ovec_test") {

  E e;
  e.y = 300.0f;
  e.z = 0.0;
  e.angle = (45 * (M_PI/180));
  double R_distance = 450.0;

  double pixel_y = 300.0;
  double pixel_z = 0.0;

  double inv_cos = 1.0 / std::cos(e.angle);
  double pow_inv_cos = inv_cos * inv_cos;

  //0vec parameters
  double vec_o[3];

  vec_o[0] = -(pixel_y + e.y - R_distance) * std::tan(e.angle) * pow_inv_cos;
  vec_o[1] = -(pixel_y + e.y + R_distance) * std::tan(e.angle) * pow_inv_cos;
  vec_o[2] = -(pixel_y + e.y) * inv_cos * (1.0 + (2.0 * std::tan(e.angle) * std::tan(e.angle)));

  CHECK(vec_o[0] == Approx(-300.0));
  CHECK(vec_o[1] == Approx(-2100.0));
  CHECK(vec_o[2] == Approx(-1800 * std::sqrt(2)));

}

TEST_CASE("avec", "avec_test") {

  E e;
  e.y = 300.0f;
  e.z = 0.0;
  e.angle = (45 * (M_PI/180));
  double R_distance = 450.0;

  double pixel_y = 300.0;
  double pixel_z = 0.0;

  double inv_cos = 1.0 / std::cos(e.angle);
  double pow_inv_cos = inv_cos * inv_cos;

  //0vec parameters
  double vec_a[3];

  vec_a[0] = -(pixel_y + e.y - R_distance) * std::tan(e.angle) * pow_inv_cos;
  vec_a[1] = -(pixel_y + e.y + R_distance) * std::tan(e.angle) * pow_inv_cos;
  vec_a[2] = -2.0 * (pixel_y + e.y) *  (inv_cos * std::tan(e.angle));

  CHECK(vec_a[0] == Approx(-300.0));
  CHECK(vec_a[1] == Approx(-2100.0));
  CHECK(vec_a[2] == Approx(-1200 * std::sqrt(2)));

}

TEST_CASE("bvec", "bvec_test") {

  E e;
  e.y = 300.0f;
  e.z = 0.0;
  e.angle = (45 * (M_PI/180));
  double R_distance = 450.0;

  double pixel_y = 300.0;
  double pixel_z = 0.0;

  double inv_cos = 1.0 / std::cos(e.angle);
  double vec_b[3];

  vec_b[0] = pixel_z - (pixel_y * std::tan(e.angle));
  vec_b[1] = pixel_z - (pixel_y * std::tan(e.angle));
  vec_b[2] = -2.0 * pixel_y * inv_cos;

  CHECK(vec_b[0] == Approx(-300.0));
  CHECK(vec_b[1] == Approx(-300.0));
  CHECK(vec_b[2] == Approx(-600 * std::sqrt(2)));

}

/*
TEST_CASE("kernel probability", "kernel_test") {

  E e;
  e.y = 300.0f;
  e.z = 0.0;
  e.angle = (45 * (M_PI/180));
  double R_distance = 450.0;
  double S_length = 450.0;
  double pixel_size = 5.0;
  int n_pixels = S_length/pixel_size;
  double sigma = 10.0;
  double dl = 63;

  spet_reconstruction<double> reconstruction(R_distance, S_length,
                                             n_pixels, pixel_size, sigma, dl);
  //reconstruction.load_input(fn);

  double y_1 = 0.0;
  double z_1 = 0.0;
  double angle = (45 * (M_PI/180));
  double z_u = -200.0f;
  double z_d = 200.0f;
  double t = reconstruction.get_event_tan(z_u, z_d);

  std::pair<int, int> p = reconstruction.pixel_center(y_1,z_1);

  std::cout << "KERNEL: " << reconstruction.kernel(y_1, z_1, angle,p)
            << std::endl;

}
*/

