#include <iostream>
#include <fstream>

#include "catch.hpp"

#include "pixel.h"
#include "lor.h"
#include "detector_ring.h"

#include "matrix_pixel_major.h"

TEST_CASE("lor/init", "LOR initalisation") {

  LOR<> lor(9, 7);

  CHECK(lor.first == 9);
  CHECK(lor.second == 7);

  CHECK(lor.index() == 9 * (9 + 1) / 2 + 7);
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
  MatrixPixelMajor<Pixel<>, LOR<>> matrix(128, 140);
}

TEST_CASE("pix_major_system_matrix/add", "add one element") {
  DetectorRing<> dr(140, 0.450, 0.006, 0.020);
  MatrixPixelMajor<Pixel<>, LOR<>> matrix(128, 140);

  LOR<> lor(9, 7);
  matrix.hit_lor(lor, 13);
  matrix.compact_pixel_index(13);

  auto hits = matrix.lor_hits_at_pixel_index(lor, 13);
  CHECK(hits == 1);
  CHECK(matrix.size() == 1);
  CHECK(matrix.n_lors_at_pixel_index(13) == 1);

  hits = matrix.lor_hits_at_pixel_index(lor, 12);
  CHECK(hits == 0);
  hits = matrix.lor_hits_at_pixel_index(LOR<>(9, 8), 13);
  CHECK(hits == 0);
}

TEST_CASE("pix_major_system_matrix/add_twice", "add one element twice") {
  DetectorRing<> dr(140, 0.450, 0.006, 0.020);
  MatrixPixelMajor<Pixel<>, LOR<>> matrix(128, 140);

  LOR<> lor(9, 7);
  matrix.hit_lor(lor, 13);
  matrix.hit_lor(lor, 13);
  matrix.compact_pixel_index(13);

  auto hits = matrix.lor_hits_at_pixel_index(lor, 13);
  CHECK(hits == 2);
  CHECK(matrix.size() == 1);
  CHECK(matrix.n_lors_at_pixel_index(13) == 1);

}

TEST_CASE("pix_major_system_matrix/add_to_all",
          "add one element to all pixels") {
  DetectorRing<> dr(140, 0.450, 0.006, 0.020);
  MatrixPixelMajor<Pixel<>, LOR<>> matrix(128, 140);

  LOR<> lor(9, 7);
  for (int i_pixel = 0; i_pixel < matrix.n_pixels(); ++i_pixel) {
    matrix.hit_lor(lor, i_pixel);
    matrix.compact_pixel_index(i_pixel);
  }

  for (int i_pixel = 0; i_pixel < matrix.n_pixels(); ++i_pixel) {
    auto hits = matrix.lor_hits_at_pixel_index(lor, i_pixel);
    CHECK(hits == 1);
    CHECK(matrix.size() == matrix.n_pixels());
    CHECK(matrix.n_lors_at_pixel_index(i_pixel) == 1);
  }

}

TEST_CASE("pix_major_system_matrix/to_sparse", "flatten") {
  DetectorRing<> dr(140, 0.450, 0.006, 0.020);
  MatrixPixelMajor<Pixel<>, LOR<>> matrix(128, 140);

  LOR<> lor(9, 7);
  for (int i_pixel = 0; i_pixel < matrix.n_pixels(); ++i_pixel) {
    matrix.hit_lor(lor, i_pixel);
    matrix.compact_pixel_index(i_pixel);
  }

  for (int i_pixel = 0; i_pixel < matrix.n_pixels(); ++i_pixel) {
    auto hits = matrix.lor_hits_at_pixel_index(lor, i_pixel);
    CHECK(hits == 1);
    CHECK(matrix.size() == matrix.n_pixels());
    CHECK(matrix.n_lors_at_pixel_index(i_pixel) == 1);
  }

  auto sparse = matrix.to_sparse();
  sparse.sort_by_lor();

  for (int i_pixel = 0; i_pixel < matrix.n_pixels(); ++i_pixel) {
    CHECK(sparse[i_pixel].lor == lor);
  }

}
