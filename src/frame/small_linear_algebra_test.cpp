#include<vector>

#include "catch.hpp"

#include"small_linear_algebra.h"

typedef float FLOAT;

TEST_CASE("Slice/01", "Slice") {
  FLOAT rep[8]={0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
  Slice<FLOAT> slice=Slice<FLOAT>(rep,0,1);
  for(int i=0;i<8;i++) 
    REQUIRE(slice[i]==Approx(rep[i]));

}


TEST_CASE("Slice/51", "Slice") {
  FLOAT rep[8]={0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
  Slice<FLOAT> slice=Slice<FLOAT>(rep,5,1);
  for(int i=0;i<3;i++) 
    REQUIRE(slice[i]==Approx(rep[i+5]));

}

TEST_CASE("Slice/32", "Slice") {
  FLOAT rep[8]={0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
  Slice<FLOAT> slice=Slice<FLOAT>(rep,3,2);
  
  REQUIRE(slice[0]==Approx(rep[3]));
  REQUIRE(slice[1]==Approx(rep[5]));
  REQUIRE(slice[2]==Approx(rep[7]));
  

}

TEST_CASE("Matrix/Square", "Square") {

  Matrix<3, 3, FLOAT> matrix;

  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      matrix[i][j]=3.0*i+j;
  Slice<FLOAT> rep =  matrix[0];
  for(int i=0;i<9;i++)
    REQUIRE(rep[i]==Approx(i));
}


TEST_CASE("Matrix/Square/Init", "Square/Init") {

  FLOAT data[9]={0.1, 0.2, 0.3,
                 0.4, 0.5, 0.6,
                 0.7, 0.8, 0.9};

  std::vector<FLOAT> vec(&data[0],&data[0]+9);

  Matrix<3, 3, FLOAT> matrix(data);

  Slice<FLOAT> rep =  matrix[0];
  for(int i=0;i<9;i++)
    REQUIRE(rep[i]==Approx(data[i]));

  Matrix<3, 3, FLOAT> matrixIt(vec.begin());

  Slice<FLOAT> repIt =  matrixIt[0];
  for(int i=0;i<9;i++)
    REQUIRE(repIt[i]==Approx(data[i]));
}


TEST_CASE("Vector/Init", "Init") {

  FLOAT data[9]={0.1, 0.2, 0.3,
                 0.4, 0.5, 0.6,
                 0.7, 0.8, 0.9};

  std::vector<FLOAT> vec(&data[0],&data[0]+9);

  Vector<9, FLOAT> vector(data);


  for(int i=0;i<9;i++)
    REQUIRE(vector[i]==Approx(data[i]));

  Vector<9, FLOAT> vectorIt(vec.begin());


  for(int i=0;i<9;i++)
    REQUIRE(vectorIt[i]==Approx(data[i]));
}

TEST_CASE("Form", "Form") {

  FLOAT data[6]={0.1, 0.2, 0.4, 0.5, 0.7, 0.8};
  Matrix<3,2,FLOAT> mat(data);
  
  FLOAT a[3]={0.7, 0.2, 0.4};
  FLOAT b[2]={0.43, 0.7};

  Vector<3,FLOAT> av(a);
  Vector<2,FLOAT> bv(b);
 
  REQUIRE(form(av,mat,bv)==Approx(0.5769));
  
}
