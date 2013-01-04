#include<iostream>
#include<fstream>

#include <catch.hpp>

#include"simple_lor.h"
#include"detector_ring.h"

#include"pix_major_matrix.h"

TEST_CASE("simple_lor/init", "simple_lor initalisation") {

  SimpleLor::init(100);
  SimpleLor lor(9,7);

  CHECK( lor.first  == 9);
  CHECK( lor.second == 7);

  CHECK( lor.n_lors()== 100*(100+1)/2);
  CHECK (SimpleLor::t_index(lor) == 9*(9+1)/2+7);
}


TEST_CASE("simple_lor/iterator", "simple_lor iterator") {

  SimpleLor::init(10);

  int count=0;
  for(auto l_it=SimpleLor::begin();l_it!=SimpleLor::end();++l_it) {
    count++;
  }
  CHECK( count== 10*(10-1)/2);
}

TEST_CASE("pix_major_system_matrix/init", "simple pix matrix test") {
  SimpleLor::init(140);
  detector_ring<> dr(140,0.450,0.006,0.020);
  PixelMajorSystemMatrix<SimpleLor> matrix(dr,128,0.005);
}


TEST_CASE("pix_major_system_matrix/add", "add one element") {
  SimpleLor::init(140);
  detector_ring<> dr(140,0.450,0.006,0.020);
  PixelMajorSystemMatrix<SimpleLor> matrix(dr,128,0.005);

  SimpleLor lor(9,7);
  matrix.add_to_t_matrix(lor,13);
  matrix.finalize_pixel(13);

  auto hit=matrix.t_get_element(lor,13);
  CHECK( hit == 1);
  CHECK( matrix.n_entries() == 1);
  CHECK( matrix.n_lors(13) == 1);

  hit=matrix.t_get_element(lor,12);
  CHECK( hit == 0);
  hit=matrix.t_get_element(SimpleLor(9,8),13);
  CHECK( hit == 0);
}


TEST_CASE("pix_major_system_matrix/add_twice", "add one element twice") {
  SimpleLor::init(140);
  detector_ring<> dr(140,0.450,0.006,0.020);
  PixelMajorSystemMatrix<SimpleLor> matrix(dr,128,0.005);

  SimpleLor lor(9,7);
  matrix.add_to_t_matrix(lor,13);
  matrix.add_to_t_matrix(lor,13);
  matrix.finalize_pixel(13);

  auto hit=matrix.t_get_element(lor,13);
  CHECK( hit == 2);
  CHECK( matrix.n_entries() == 1);
  CHECK( matrix.n_lors(13) == 1);


}



TEST_CASE("pix_major_system_matrix/add_to_all", "add one element to all pixels") {
  SimpleLor::init(140);
  detector_ring<> dr(140,0.450,0.006,0.020);
  PixelMajorSystemMatrix<SimpleLor> matrix(dr,128,0.005);

    SimpleLor lor(9,7);
    for(int p=0;p<matrix.total_n_pixels();++p) {
        matrix.add_to_t_matrix(lor,p);
        matrix.finalize_pixel(p);
  }

    for(int p=0;p<matrix.total_n_pixels();++p) {
      auto hit=matrix.t_get_element(lor,p);
      CHECK( hit == 1);
      CHECK( matrix.n_entries() == matrix.total_n_pixels());
      CHECK( matrix.n_lors(p) == 1);
  }

}
