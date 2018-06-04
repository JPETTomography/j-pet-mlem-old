#include "util/test.h"

#include "3d/geometry/matrix.h"
#include "3d/geometry/vector.h"

#include "common/types.h"

using Matrix = PET3D::Matrix<F>;
using Vector = PET3D::Vector<Matrix::F>;

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
    REQUIRE(numbers[i] == Approx(i + 1).epsilon(1e-7));
}

TEST("3d/geometry/matrix/vector_multiplication") {
  Matrix mat{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
  Vector vec = { 1, 2, 3 };
  CHECK(vec.x == 1.0_e7);
  CHECK(vec.y == 2.0_e7);
  CHECK(vec.z == 3.0_e7);

  Vector res = mat * vec;

  REQUIRE(res.x == 14.0_e7);
  REQUIRE(res.y == 32.0_e7);
  REQUIRE(res.z == 50.0_e7);
}

TEST("3d/geometry/matrix/arithmetic_assignment_operators") {
  Matrix rhs{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
  {
    Matrix mat{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    mat += rhs;
    for (int i = 0; i < 9; i++) {
      REQUIRE(mat[i] == Approx(2 * (i + 1)).epsilon(1e-7));
    }
  }
  {
    Matrix mat{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    mat -= rhs;
    for (int i = 0; i < 9; i++) {
      REQUIRE(mat[i] == Approx(0).epsilon(1e-7));
    }
  }
}

TEST("3d/geometry/matrix/arithmetic_operators") {
  Matrix rhs{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
  Matrix lhs{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
  {
    Matrix mat;

    mat = lhs + rhs;
    for (int i = 0; i < 9; i++) {
      REQUIRE(mat[i] == Approx(2 * (i + 1)).epsilon(1e-7));
    }
  }
  {
    Matrix mat;
    mat = lhs - rhs;
    for (int i = 0; i < 9; i++) {
      REQUIRE(mat[i] == Approx(0).epsilon(1e-7));
    }
  }
}

TEST("3d/geometry/matrix/scalar_multiplication") {
  Matrix rhs{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
  {
    Matrix mat;

    mat = rhs * F(3);
    for (int i = 0; i < 9; i++) {
      REQUIRE(mat[i] == Approx(3 * (i + 1)).epsilon(1e-7));
    }
  }
  {
    Matrix mat;
    mat = F(3) * rhs;
    for (int i = 0; i < 9; i++) {
      REQUIRE(mat[i] == Approx(3 * (i + 1)).epsilon(1e-7));
    }
  }
}

TEST("3d/geometry/matrix/rotation") {
  {
    auto rotated = Matrix::rotate(Vector(1, 0, 0), 90._deg) * Vector(0, 1, 0);
    REQUIRE((rotated - Vector(0, 0, 1)).length() == 0._e7);
  }
  {
    auto rotated = Matrix::rotate(Vector(1, 0, 0), 180._deg) * Vector(0, 1, 0);
    REQUIRE((rotated - Vector(0, -1, 0)).length() == 0._e7);
  }
  {
    auto rotated = Matrix::rotate(Vector(0, 0, 1), 90._deg) * Vector(0, 1, 0);
    REQUIRE((rotated - Vector(-1, 0, 0)).length() == 0._e7);
  }
  {
    auto rotated = Matrix::rotate(Vector(0, 0, 1), 180._deg) * Vector(0, 1, 0);
    REQUIRE((rotated - Vector(0, -1, 0)).length() == 0._e7);
  }
}
