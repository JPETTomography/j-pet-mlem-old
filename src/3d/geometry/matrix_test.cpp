
#include "util/test.h"

#include "3d/geometry/matrix.h"
#include "3d/geometry/vector.h"

using Matrix = PET3D::Matrix<float>;

TEST("3d/geometry/matrix/initialisation") {

  Matrix mat;

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      REQUIRE(mat(i, j) == 0.0_e7);

  Matrix one = Matrix::identity();

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      if (i == j)
        REQUIRE(one(i, j) == 1.0_e7);
      else
        REQUIRE(one(i, j) == 0.0_e7);
    }

  Matrix numbers{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
  for (int i = 0; i < 9; i++)
    REQUIRE(numbers(i) == Approx(i + 1).epsilon(1e-7));
}

TEST("3d/geometry/matrix/vector_multiplication") {
  using Vector = PET3D::Vector<Matrix::F>;

  Matrix mat{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
  Vector vec = { 1, 2, 3 };
  CHECK(vec.x == 1.0_e7);
  CHECK(vec.y == 2.0_e7);
  CHECK(vec.z == 3.0_e7);

  Vector res = mat * vec;

  REQUIRE(res.x==14.0_e7);
  REQUIRE(res.y==32.0_e7);
  REQUIRE(res.z==50.0_e7);

}