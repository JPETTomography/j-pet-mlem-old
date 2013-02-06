#include <iostream>
#include <fstream>

#include "catch.hpp"

#include "lor.h"
#include "detector_ring.h"

#include "matrix_pixel_major.h"

TEST_CASE("lor/init", "LOR initalisation") {

  LOR<> lor(9, 7);

  CHECK(lor.first == 9);
  CHECK(lor.second == 7);

  CHECK(lor.t_index() == 9 * (9 + 1) / 2 + 7);
}

TEST_CASE("lor/iterator", "LOR iterator") {

  int count = 0;
  for (auto lor = LOR<>(); lor != LOR<>::end_for_detectors(10); ++lor) {
    count++;
  }
  CHECK(count == 10 * (10 + 1) / 2);
}

TEST_CASE("pix_major_system_matrix/init", "simple pix matrix test") {
  DetectorRing<> dr(140, 0.450, 0.006, 0.020);
  MatrixPixelMajor<LOR<>> matrix(128, 140);
}

TEST_CASE("pix_major_system_matrix/add", "add one element") {
  DetectorRing<> dr(140, 0.450, 0.006, 0.020);
  MatrixPixelMajor<LOR<>> matrix(128, 140);

  LOR<> lor(9, 7);
  matrix.add_to_t_matrix(lor, 13);
  matrix.compact_pixel_index(13);

  auto hit = matrix.t_get_element(lor, 13);
  CHECK(hit == 1);
  CHECK(matrix.size() == 1);
  CHECK(matrix.n_lors(13) == 1);

  hit = matrix.t_get_element(lor, 12);
  CHECK(hit == 0);
  hit = matrix.t_get_element(LOR<>(9, 8), 13);
  CHECK(hit == 0);
}

TEST_CASE("pix_major_system_matrix/add_twice", "add one element twice") {
  DetectorRing<> dr(140, 0.450, 0.006, 0.020);
  MatrixPixelMajor<LOR<>> matrix(128, 140);

  LOR<> lor(9, 7);
  matrix.add_to_t_matrix(lor, 13);
  matrix.add_to_t_matrix(lor, 13);
  matrix.compact_pixel_index(13);

  auto hit = matrix.t_get_element(lor, 13);
  CHECK(hit == 2);
  CHECK(matrix.size() == 1);
  CHECK(matrix.n_lors(13) == 1);

}

TEST_CASE("pix_major_system_matrix/add_to_all",
          "add one element to all pixels") {
  DetectorRing<> dr(140, 0.450, 0.006, 0.020);
  MatrixPixelMajor<LOR<>> matrix(128, 140);

  LOR<> lor(9, 7);
  for (int p = 0; p < matrix.total_n_pixels(); ++p) {
    matrix.add_to_t_matrix(lor, p);
    matrix.compact_pixel_index(p);
  }

  for (int p = 0; p < matrix.total_n_pixels(); ++p) {
    auto hit = matrix.t_get_element(lor, p);
    CHECK(hit == 1);
    CHECK(matrix.size() == matrix.total_n_pixels());
    CHECK(matrix.n_lors(p) == 1);
  }

}

TEST_CASE("pix_major_system_matrix/to_pairs", "flatten") {
  DetectorRing<> dr(140, 0.450, 0.006, 0.020);
  MatrixPixelMajor<LOR<>> matrix(128, 140);

  LOR<> lor(9, 7);
  for (int p = 0; p < matrix.total_n_pixels(); ++p) {
    matrix.add_to_t_matrix(lor, p);
    matrix.compact_pixel_index(p);
  }

  for (int p = 0; p < matrix.total_n_pixels(); ++p) {
    auto hit = matrix.t_get_element(lor, p);
    CHECK(hit == 1);
    CHECK(matrix.size() == matrix.total_n_pixels());
    CHECK(matrix.n_lors(p) == 1);
  }

  matrix.make_list();
  matrix.sort_pairs_by_lors();

  for (int p = 0; p < matrix.total_n_pixels(); ++p) {
    CHECK(std::get<0>(matrix[p]) == lor);
  }

}
