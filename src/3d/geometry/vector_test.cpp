#include "util/test.h"

#include "3d/geometry/vector.h"

#include "common/types.h"

using Vector = PET3D::Vector<F>;

TEST("3d/geometry/vector/init") {

  Vector vec(1, 2, 3);

  CHECK(vec.x == 1.0_e7);
  CHECK(vec.y == 2.0_e7);
  CHECK(vec.z == 3.0_e7);
}

TEST("3d/geometry/vector/arithmetic_assignement") {
  {
    Vector vec1(1, 2, 3);
    Vector vec2(0.1, 0.2, 0.3);
    vec1 += vec2;

    CHECK(vec1.x == 1.1_e7);
    CHECK(vec1.y == 2.2_e7);
    CHECK(vec1.z == 3.3_e7);
  }
  {
    Vector vec1(1, 2, 3);
    Vector vec2(0.1, 0.2, 0.3);
    vec1 -= vec2;

    CHECK(vec1.x == 0.9_e7);
    CHECK(vec1.y == 1.8_e7);
    CHECK(vec1.z == 2.7_e7);
  }
  {
    Vector vec1(1, 2, 3);
    vec1 *= 2;

    CHECK(vec1.x == 2.0_e7);
    CHECK(vec1.y == 4.0_e7);
    CHECK(vec1.z == 6.0_e7);
  }
  {
    Vector vec1(1, 2, 3);
    vec1 /= 2;

    CHECK(vec1.x == 0.5_e7);
    CHECK(vec1.y == 1.0_e7);
    CHECK(vec1.z == 1.5_e7);
  }
}

TEST("3d/geometry/vector/arithmetics") {
  {
    Vector lhs(1, 2, 3);
    Vector rhs(0.1, 0.2, 0.3);
    Vector vec = lhs + rhs;

    CHECK(vec.x == 1.1_e7);
    CHECK(vec.y == 2.2_e7);
    CHECK(vec.z == 3.3_e7);
  }
  {
    Vector lhs(1, 2, 3);
    Vector rhs(0.1, 0.2, 0.3);
    Vector vec = lhs - rhs;

    CHECK(vec.x == 0.9_e7);
    CHECK(vec.y == 1.8_e7);
    CHECK(vec.z == 2.7_e7);
  }
}

TEST("3d/geometry/vector/logical") {
  {
    Vector lhs(1, 2, 3);
    Vector rhs(1, 2, 3);

    CHECK(lhs == rhs);
    CHECK(!(lhs != rhs));
  }
  {
    Vector lhs(1, 2, 4);
    Vector rhs(1, 2, 3);

    CHECK(!(lhs == rhs));
    CHECK((lhs != rhs));
  }
  {
    Vector lhs(1, 5, 3);
    Vector rhs(1, 2, 3);

    CHECK(!(lhs == rhs));
    CHECK((lhs != rhs));
  }
  {
    Vector lhs(5, 2, 3);
    Vector rhs(1, 2, 3);

    CHECK(!(lhs == rhs));
    CHECK((lhs != rhs));
  }
}
